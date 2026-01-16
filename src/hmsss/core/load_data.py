from __future__ import annotations

import sys
import hashlib
import shutil
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen, Request


class DataBootstrapError(RuntimeError):
    pass


def _project_root_from_this_file(this_file: Path, *, marker_dir: str = "src") -> Path:
    """
    Robust: .../HMSS2/src/... -> Root ist .../HMSS2
    """
    p = this_file.resolve()
    for parent in [p] + list(p.parents):
        if parent.name == marker_dir:
            return parent.parent
    raise DataBootstrapError(f"Could not infer project root from: {this_file}")


def _download_file(
    url: str,
    out_path: Path,
    *,
    user_agent: str = "HMSS2/1.0",
    chunk_size: int = 1024 * 1024,  # 1 MB
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    req = Request(url, headers={"User-Agent": user_agent})
    with urlopen(req) as r, open(out_path, "wb") as f:
        total_size = r.headers.get("Content-Length")
        total_size = int(total_size) if total_size is not None else None

        downloaded = 0

        for chunk in iter(lambda: r.read(chunk_size), b""):
            f.write(chunk)
            downloaded += len(chunk)

            if total_size:
                pct = (downloaded / total_size) * 100
                sys.stdout.write(
                    f"\r[Download] {downloaded / 1e6:7.1f} / {total_size / 1e6:7.1f} MB "
                    f"({pct:3.0f}%)"
                )
            else:
                sys.stdout.write(f"\r[Download] {downloaded / 1e6:7.1f} MB")

            sys.stdout.flush()

    sys.stdout.write("\n")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dir_with_marker_from_zip(
    *,
    zip_url: str,
    expected_sha256: str | None = None,
    project_root: Path | None = None,
    target_dirname: str,
    marker_name: str = ".ok",
    user_agent: str = "HMSS2/1.0",
) -> Path:
    """
    Stellt sicher, dass <project_root>/<target_dirname> vollständig initialisiert ist.

    "Vollständig" = target_dir existiert UND Marker-Datei existiert.

    ZIP-Layouts unterstützt:
      A) ZIP enthält top-level '<target_dirname>/...'
      B) ZIP enthält direkt den Inhalt von '<target_dirname>/' (z.B. 'HMMs/...')

    Verhalten wenn target_dir existiert, aber Marker fehlt:
      - target_dir wird als unvollständig betrachtet
      - target_dir wird nach <target_dirname>.__backup__<timestamp> verschoben
      - neues target_dir wird aus dem ZIP erstellt
      - Marker wird am Ende angelegt
    """
    if project_root is None:
        project_root = _project_root_from_this_file(Path(__file__))

    target_dir = project_root / target_dirname
    marker = target_dir / marker_name

    # Skip nur wenn Marker existiert
    if target_dir.is_dir() and marker.is_file():
        return target_dir

    if not target_dir.exists():
        sys.stdout.write(
            f"[SETUP] Archive is not initialized: {target_dir}\n"
            f"[SETUP] Downloading from: {zip_url}\n"
        )
    else:
        sys.stdout.write(
            f"[SETUP] Archive is incomplete: {target_dir}\n"
            f"[SETUP] Downloading from: {zip_url}\n"
        )

    project_root.mkdir(parents=True, exist_ok=True)

    # Backup, falls unvollständig vorhanden
    backup_dir: Path | None = None
    if target_dir.exists():
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = project_root / f"{target_dirname}.__backup__{ts}"
        i = 1
        while backup_dir.exists():
            backup_dir = project_root / f"{target_dirname}.__backup__{ts}_{i}"
            i += 1
        target_dir.replace(backup_dir)

    tmp_root = Path(tempfile.mkdtemp(prefix=f"hmss2_{target_dirname}_bootstrap_"))
    try:
        zip_path = tmp_root / f"{target_dirname}.zip"
        extract_dir = tmp_root / "extract"
        extract_dir.mkdir(parents=True, exist_ok=True)

        # Download
        _download_file(zip_url, zip_path, user_agent=user_agent)

        # Optional integrity check
        if expected_sha256 is not None:
            got = _sha256(zip_path)
            if got.lower() != expected_sha256.lower():
                raise DataBootstrapError(
                    f"SHA256 mismatch for {zip_path}: expected {expected_sha256}, got {got}"
                )

        # Extract
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_dir)

        # Quelle bestimmen (Layout A oder B)
        candidate_a = extract_dir / target_dirname
        source_dir = candidate_a if candidate_a.is_dir() else extract_dir

        # Staging copy -> atomarer replace
        staging = project_root / f".{target_dirname}_staging_tmp"
        if staging.exists():
            shutil.rmtree(staging)
        shutil.copytree(source_dir, staging)

        staging.replace(target_dir)

        if not target_dir.is_dir():
            raise DataBootstrapError(
                f"Target directory was not created after extraction: {target_dir}"
            )

        # Marker am Ende
        marker.touch()
        return target_dir

    except zipfile.BadZipFile as e:
        raise DataBootstrapError(f"Downloaded file is not a valid zip: {e}") from e
    except Exception as e:
        # Backup zurückholen, falls möglich
        if (
            (not target_dir.exists())
            and (backup_dir is not None)
            and backup_dir.exists()
        ):
            try:
                backup_dir.replace(target_dir)
            except Exception:
                pass
        raise DataBootstrapError(f"Failed to bootstrap {target_dirname}: {e}") from e
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)


def ensure_data_dir_with_marker(
    *,
    data_zip_url: str,
    expected_sha256: str | None = None,
    project_root: Path | None = None,
    marker_name: str = ".data_ok",
) -> Path:
    return ensure_dir_with_marker_from_zip(
        zip_url=data_zip_url,
        expected_sha256=expected_sha256,
        project_root=project_root,
        target_dirname="data",
        marker_name=marker_name,
    )


def ensure_gpkg_dir_with_marker(
    *,
    gpkg_zip_url: str,
    expected_sha256: str | None = None,
    project_root: Path | None = None,
    marker_name: str = ".gpkg_ok",
) -> Path:
    return ensure_dir_with_marker_from_zip(
        zip_url=gpkg_zip_url,
        expected_sha256=expected_sha256,
        project_root=project_root,
        target_dirname="gpkg",
        marker_name=marker_name,
    )
