from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import zipfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from urllib.request import Request, urlopen

PACKAGE_MANIFEST = Path(__file__).with_name("data_manifest.json")
INSTALLED_MANIFEST = ".unisul_data_manifest.json"
CONFIG_FILENAME = "config.json"


class DataManagementError(RuntimeError):
    """Raised when UNISUL reference-data configuration or installation fails."""


def user_config_dir() -> Path:
    if os.name == "nt":
        return Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support"
    return Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))


def user_data_dir() -> Path:
    if os.name == "nt":
        return Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support"
    return Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))


def config_file() -> Path:
    return user_config_dir() / "unisul" / CONFIG_FILENAME


def _configured_data_path() -> Path | None:
    path = config_file()
    if not path.is_file():
        return None
    try:
        contents = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DataManagementError(
            f"Cannot read UNISUL configuration {path}: {exc}"
        ) from exc
    value = contents.get("data_dir") if isinstance(contents, dict) else None
    if not value:
        return None
    configured = Path(value).expanduser()
    if not configured.is_absolute():
        raise DataManagementError(
            f"Configured data_dir must be an absolute path: {path}"
        )
    return configured


def _resolve_data_root(cli_path: str | os.PathLike | None = None) -> tuple[Path, str]:
    """Resolve the directory that directly contains the complete UNISUL data tree."""
    if cli_path is not None:
        return Path(cli_path).expanduser().resolve(), "CLI"
    env = os.environ.get("UNISUL_DATA")
    if env:
        return Path(env).expanduser().resolve(), "Environment"
    configured = _configured_data_path()
    if configured is not None:
        return configured.resolve(), "Config"
    return (user_data_dir() / "unisul" / "data").resolve(), "Default"


def resolve_data_root(cli_path: str | os.PathLike | None = None) -> Path:
    """Resolve data path in CLI > UNISUL_DATA > saved config > default order."""
    return _resolve_data_root(cli_path)[0]


def data_root_source(cli_path: str | os.PathLike | None = None) -> str:
    return _resolve_data_root(cli_path)[1]


def load_package_manifest() -> dict:
    try:
        manifest = json.loads(PACKAGE_MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DataManagementError(f"Cannot read packaged data manifest: {exc}") from exc
    if not isinstance(manifest, dict) or "archive" not in manifest:
        raise DataManagementError(f"Invalid packaged data manifest: {PACKAGE_MANIFEST}")
    return manifest


def _installed_manifest(root: Path) -> dict:
    path = root / INSTALLED_MANIFEST
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DataManagementError(
            f"Cannot read installed data manifest {path}: {exc}"
        ) from exc
    return value if isinstance(value, dict) else {}


def required_paths() -> list[str]:
    return list(load_package_manifest().get("required_paths", []))


def missing_required_paths(root: Path) -> list[str]:
    return [name for name in required_paths() if not (root / name).exists()]


def data_is_present(root: Path) -> bool:
    return root.is_dir() and not missing_required_paths(root)


def require_pipeline_data(config) -> None:
    root = Path(config.paths.data_root or config.paths.data).expanduser().resolve()
    missing = missing_required_paths(root)
    if missing:
        details = "\n".join(f"  - {name}" for name in missing)
        raise DataManagementError(
            "Required UNISUL reference data are missing.\n\n"
            f"Resolved data directory:\n  {root}\n\n"
            f"Missing required entries:\n{details}\n\n"
            "Install the complete reference-data archive with:\n"
            "  unisul data download\n\n"
            "or select an existing complete data directory using:\n"
            "  --data-dir PATH\n"
            "  export UNISUL_DATA=PATH\n"
            "  unisul data set-path PATH"
        )


def set_data_path(path: str | os.PathLike) -> Path:
    target = Path(path).expanduser().resolve()
    cfg = config_file()
    cfg.parent.mkdir(parents=True, exist_ok=True)
    temporary = cfg.with_name(f".{cfg.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps({"data_dir": str(target)}, indent=2) + "\n", encoding="utf-8"
    )
    os.replace(temporary, cfg)
    return cfg


@contextmanager
def _download_lock(root: Path):
    parent = root.parent
    parent.mkdir(parents=True, exist_ok=True)
    lock = parent / f".{root.name}.unisul-download.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise DataManagementError(
            f"A UNISUL data download is already running (lock: {lock}). "
            "Remove the lock only after confirming no download process is active."
        ) from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid={os.getpid()}\n")
        yield
    finally:
        lock.unlink(missing_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_archive_path(staging: Path, member_name: str) -> Path:
    normalized = PurePosixPath(member_name)
    if normalized.is_absolute() or ".." in normalized.parts:
        raise DataManagementError(f"Unsafe archive path: {member_name}")
    target = (staging / Path(*normalized.parts)).resolve()
    if not target.is_relative_to(staging.resolve()):
        raise DataManagementError(f"Unsafe archive path: {member_name}")
    return target


def _extract_archive(archive: Path, staging: Path, archive_format: str) -> None:
    if archive_format == "zip":
        with zipfile.ZipFile(archive) as source:
            for member in source.infolist():
                target = _safe_archive_path(staging, member.filename)
                mode = member.external_attr >> 16
                if mode & 0o170000 == 0o120000:
                    raise DataManagementError(
                        "Symbolic links are not allowed in data archives"
                    )
                if member.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with source.open(member) as inp, target.open("wb") as out:
                        shutil.copyfileobj(inp, out)
        return

    modes = {
        "tar": "r:",
        "tar.gz": "r:gz",
        "tgz": "r:gz",
        "tar.xz": "r:xz",
        "tar.bz2": "r:bz2",
    }
    mode = modes.get(archive_format)
    if mode is None:
        raise DataManagementError(f"Unsupported archive format: {archive_format}")

    with tarfile.open(archive, mode) as source:
        for member in source.getmembers():
            target = _safe_archive_path(staging, member.name)
            if not (member.isfile() or member.isdir()):
                raise DataManagementError(
                    f"Unsupported archive entry type: {member.name}"
                )
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            extracted = source.extractfile(member)
            if extracted is None:
                raise DataManagementError(
                    f"Cannot extract archive entry: {member.name}"
                )
            with extracted, target.open("wb") as out:
                shutil.copyfileobj(extracted, out)


def _select_extracted_data_root(extract_dir: Path) -> Path:
    """Accept archives containing data/... or the direct contents of data/."""
    direct = extract_dir / "data"
    if direct.is_dir() and not missing_required_paths(direct):
        return direct

    if not missing_required_paths(extract_dir):
        return extract_dir

    top_dirs = [p for p in extract_dir.iterdir() if p.is_dir()]
    top_files = [p for p in extract_dir.iterdir() if p.is_file()]
    if len(top_dirs) == 1 and not top_files and not missing_required_paths(top_dirs[0]):
        return top_dirs[0]

    missing = missing_required_paths(extract_dir)
    raise DataManagementError(
        "Downloaded archive does not contain a valid complete UNISUL data directory. "
        f"Missing expected entries include: {', '.join(missing[:8])}"
    )


def _write_installed_manifest(root: Path, package_manifest: dict) -> None:
    payload = {
        "data_version": package_manifest.get("data_version"),
        "minimum_unisul_version": package_manifest.get("minimum_unisul_version"),
        "archive_filename": package_manifest.get("archive", {}).get("filename"),
    }
    target = root / INSTALLED_MANIFEST
    temporary = root / f".{INSTALLED_MANIFEST}.{os.getpid()}.tmp"
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, target)


def download_data(root: Path, *, force: bool = False) -> Path:
    """Download and atomically install one complete UNISUL data archive."""
    root = Path(root).expanduser().resolve()
    manifest = load_package_manifest()
    archive_cfg = manifest["archive"]
    url = archive_cfg.get("url")
    checksum = archive_cfg.get("sha256")

    if not url:
        raise DataManagementError(
            "No reference-data download URL is configured in data_manifest.json."
        )

    if data_is_present(root) and not force:
        return root
    if root.exists() and not force:
        raise DataManagementError(
            f"Data directory exists but is incomplete: {root}\n"
            "Move it aside, select another --data-dir, or rerun with "
            "'unisul data download --force'."
        )

    with _download_lock(root):
        if data_is_present(root) and not force:
            return root

        root.parent.mkdir(parents=True, exist_ok=True)
        temporary_root = Path(
            tempfile.mkdtemp(prefix=".unisul-data-download.", dir=root.parent)
        )
        archive = temporary_root / archive_cfg.get("filename", "unisul-data.archive")
        partial = archive.with_name(archive.name + ".part")
        extract_dir = temporary_root / "extract"
        extract_dir.mkdir()
        backup = None

        try:
            request = Request(url, headers={"User-Agent": "UNISUL"})
            with urlopen(request, timeout=60) as response, partial.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
                output.flush()
                os.fsync(output.fileno())

            if checksum:
                actual = _sha256(partial)
                if actual.lower() != str(checksum).lower():
                    raise DataManagementError(
                        f"SHA256 mismatch: expected {checksum}, got {actual}"
                    )

            os.replace(partial, archive)
            _extract_archive(archive, extract_dir, archive_cfg["archive_format"])
            source = _select_extracted_data_root(extract_dir)

            staging = root.parent / f".{root.name}.staging.{os.getpid()}"
            if staging.exists():
                shutil.rmtree(staging)
            shutil.copytree(source, staging)
            staging.chmod(stat.S_IMODE(root.parent.stat().st_mode))

            if root.exists():
                backup = root.parent / f".{root.name}.backup.{os.getpid()}"
                if backup.exists():
                    shutil.rmtree(backup)
                os.replace(root, backup)

            os.replace(staging, root)
            _write_installed_manifest(root, manifest)

            if backup is not None:
                shutil.rmtree(backup, ignore_errors=True)
            return root
        except Exception as exc:
            if backup is not None and backup.exists() and not root.exists():
                os.replace(backup, root)
            if isinstance(exc, DataManagementError):
                raise
            raise DataManagementError(
                f"Failed to download UNISUL reference data: {exc}"
            ) from exc
        finally:
            shutil.rmtree(temporary_root, ignore_errors=True)


def _print_status(root: Path, cli_path: str | None) -> None:
    installed = _installed_manifest(root)
    missing = missing_required_paths(root)
    print(f"Data directory: {root}")
    print(f"Selected by: {data_root_source(cli_path)}")
    print(f"Directory exists: {'yes' if root.is_dir() else 'no'}")
    print(f"Complete: {'yes' if not missing else 'no'}")
    print(f"Data version: {installed.get('data_version') or 'unknown'}")
    if missing:
        print("Missing:")
        for name in missing:
            print(f"  {name}")


def data_cli(argv: list[str], cli_path: str | None = None) -> None:
    parser = argparse.ArgumentParser(prog="unisul data")
    subparsers = parser.add_subparsers(dest="command", required=True)

    download = subparsers.add_parser(
        "download", help="Download the complete UNISUL reference-data archive"
    )
    download.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing incomplete or old data directory",
    )
    subparsers.add_parser(
        "status", help="Show selected reference-data directory and status"
    )
    set_path = subparsers.add_parser(
        "set-path", help="Persist a reference-data directory"
    )
    set_path.add_argument("path")

    args = parser.parse_args(argv)
    root = resolve_data_root(cli_path)

    try:
        if args.command == "status":
            _print_status(root, cli_path)
        elif args.command == "set-path":
            cfg = set_data_path(args.path)
            print(f"Saved reference-data path to {cfg}")
            print(f"Data directory: {Path(args.path).expanduser().resolve()}")
        elif args.command == "download":
            target = download_data(root, force=args.force)
            print(f"Installed complete UNISUL reference data in {target}")
    except DataManagementError as exc:
        parser.exit(2, f"unisul data: error: {exc}\n")
