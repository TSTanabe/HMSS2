# tests/test_config_facade.py
from pathlib import Path
import pytest

from hmsss.core.config import Config, PathsCfg  # basiert auf deiner config.py

@pytest.fixture
def cfg(tmp_path) -> Config:
    """
    Baut eine minimale HMSS2-Verzeichnisstruktur im Temp-Verzeichnis und
    erzeugt daraus ein valides Config-Objekt.
    """
    root = tmp_path / "HMSS2"
    (root / "bin").mkdir(parents=True)
    (root / "data" / "HMMlib").mkdir(parents=True)
    (root / "data" / "RefSeqs").mkdir(parents=True)
    (root / "results").mkdir(parents=True)
    (root / "src" / "hmsss").mkdir(parents=True)

    paths = PathsCfg(
        root=str(root),
        bin=str(root / "bin"),
        data=str(root / "data"),
        hmms=str(root / "data" / "HMMlib"),
        refseq=str(root / "data" / "RefSeqs"),
        results=str(root / "results"),
        package=str(root / "src" / "hmsss"),
    )
    c = Config(paths=paths)
    # validiert u. a. die Basisverzeichnisse
    c.validate()
    return c


def test_facade_stage_and_lists(cfg: Config):
    # Default kommt aus CliSearchParams.stage == 0
    assert cfg.stage == 0
    # Schreiben via Fassade -> landet in cfg.cli_params.stage
    cfg.stage = 101
    assert cfg.cli_params.stage == 101

    # hmm_sets (Liste) lesen/schreiben über Fassade
    assert cfg.hmm_sets == []
    cfg.hmm_sets = ["grpA", "grpB"]
    assert cfg.cli_resources.HMM_sets == ["grpA", "grpB"]

    # Auch Mutation der Liste wird durchgereicht
    cfg.hmm_sets.append("grpC")
    assert cfg.cli_resources.HMM_sets == ["grpA", "grpB", "grpC"]

    # Weitere Facade-Beispiele (Operators)
    cfg.fetch_proteins = ["PF00001", "PF00002"]
    assert cfg.cli_ops.fetch_proteins == ["PF00001", "PF00002"]


def test_facade_result_dir_priority(cfg: Config, tmp_path):
    # Lesen: result_dir priorisiert project.result_files_directory,
    # sonst CLI-Wunsch, sonst paths.results
    default_results = Path(cfg.paths.results)
    assert Path(cfg.result_dir) == default_results

    # Schreiben: setter setzt CLI-Wunsch (cli_input.result_files_directory)
    user_dir = tmp_path / "custom_results"
    cfg.result_dir = str(user_dir)
    assert cfg.cli_input.result_files_directory == str(user_dir)
    assert Path(cfg.result_dir) == user_dir

    # Wenn das Projekt später ein finales results_dir setzt, hat das Priorität
    project_dir = tmp_path / "final_project"
    cfg.project.result_files_directory = str(project_dir)
    assert Path(cfg.result_dir) == project_dir
