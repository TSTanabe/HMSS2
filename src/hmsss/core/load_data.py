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


def ensure_data_dir_with_marker(
    *,
    data_zip_url: str,
    expected_sha256: str | None = None,
    project_root: Path | None = None,
    marker_name: str = ".data_ok",
) -> Path:
    """
    Stellt sicher, dass <project_root>/data vollständig initialisiert ist.
    "Vollständig" = data/ existiert UND Marker-Datei existiert.

    ZIP-Layouts unterstützt:
      A) ZIP enthält top-level 'data/...'
      B) ZIP enthält direkt den Inhalt von 'data/' (z.B. 'HMMs/...', 'RefSeqs/...')

    Verhalten wenn data/ existiert, aber Marker fehlt:
      - data/ wird als unvollständig betrachtet
      - data/ wird nach data.__backup__<timestamp> verschoben
      - neues data/ wird aus dem ZIP erstellt
      - Marker wird am Ende angelegt
      - Backup bleibt erhalten (damit nichts verloren geht)
    """
    if project_root is None:
        project_root = _project_root_from_this_file(Path(__file__))

    data_dir = project_root / "data"
    marker = data_dir / marker_name

    # Nur dann skippen, wenn Marker existiert
    if data_dir.is_dir() and marker.is_file():
        return data_dir

    if not data_dir.exists():
        sys.stdout.write(
            f"[SETUP] Data archive is not initialized: {data_dir} "
            f"[SETUP] Downloading data from archive.\n"
        )
    else:
        sys.stdout.write(
            f"[SETUP] Data archive is incomplete: {data_dir} "
            f"[SETUP] Downloading data from archive.\n"
        )

    project_root.mkdir(parents=True, exist_ok=True)

    # Falls data/ existiert aber unvollständig: sichern (nicht löschen)
    backup_dir: Path | None = None
    if data_dir.exists():
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = project_root / f"data.__backup__{ts}"
        # Wenn es irgendwie schon existiert (sehr selten): suffix
        i = 1
        while backup_dir.exists():
            backup_dir = project_root / f"data.__backup__{ts}_{i}"
            i += 1
        data_dir.replace(backup_dir)

    tmp_root = Path(tempfile.mkdtemp(prefix="hmss2_data_bootstrap_"))
    try:
        zip_path = tmp_root / "data.zip"
        extract_dir = tmp_root / "extract"
        extract_dir.mkdir(parents=True, exist_ok=True)

        # Download
        _download_file(data_zip_url, zip_path)

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

        # Quelle bestimmen
        candidate_a = extract_dir / "data"
        if candidate_a.is_dir():
            source_data = candidate_a
        else:
            source_data = extract_dir

        # In staging kopieren (damit replace sauber/atomar bleibt)
        staging = project_root / ".data_staging_tmp"
        if staging.exists():
            shutil.rmtree(staging)
        shutil.copytree(source_data, staging)

        # Jetzt staging -> data (data existiert hier nicht mehr, weil ggf. vorher gebackupt)
        staging.replace(data_dir)

        # Minimal sanity check
        if not data_dir.is_dir():
            raise DataBootstrapError("Data directory was not created after extraction.")

        # ✅ Marker am Ende schreiben: nur wenn alles erfolgreich war
        marker.touch()

        return data_dir

    except zipfile.BadZipFile as e:
        raise DataBootstrapError(f"Downloaded file is not a valid zip: {e}") from e
    except Exception as e:
        # Wenn etwas schiefgeht: versuche Backup zurückzuholen
        if (not data_dir.exists()) and (backup_dir is not None) and backup_dir.exists():
            try:
                backup_dir.replace(data_dir)
            except Exception:
                pass
        raise DataBootstrapError(f"Failed to bootstrap data directory: {e}") from e
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
