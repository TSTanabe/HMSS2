import os
from pathlib import Path

import pytest

# Passe die Importe ggf. an deine Paketstruktur an:
from hmsss.core.options import Hmsss
from hmsss.db.project import prepare_result_space


def test_prepare_result_space_new_project_sets_paths_and_subdirs(tmp_path: Path):
    # Arrange
    opts = Hmsss(
        location=str(tmp_path),
        result_files_directory=str(tmp_path / "results"),
        new_project=True,  # erzwinge Neuanlage unter <location>/results
        project_name="unittest",
    )

    # Act
    prepare_result_space(opts, project=opts.project_name)

    # Assert: result_files_directory wurde erzeugt und hat den Projektnamen als Suffix
    rdir = Path(opts.result_files_directory)
    assert rdir.parent == tmp_path / "results"
    assert rdir.name.endswith("_unittest")

    # Verzeichnisse existieren physisch
    expected_subdirs = {
        "Hit_list",
        "Sequences",
        "Filtered_hits",
        "Collinear_syntenic_blocks",
    }
    for sub in expected_subdirs:
        p = rdir / sub
        assert p.is_dir(), f"Expected subdir missing: {p}"

    # options-Felder korrekt gesetzt
    assert opts.database_directory == str(rdir / "database.db")
    assert opts.fasta_initial_hit_directory == str(rdir / "Hit_list")
    assert opts.fasta_output_directory == str(rdir / "Sequences")
    assert opts.cross_check_directory == str(rdir / "Filtered_hits")
    assert opts.Csb_directory == str(rdir / "Collinear_syntenic_blocks")

    # glob_report wird gesetzt, wenn vorher None/fehlend
    assert opts.glob_report == str(rdir / "global_report.cat_hmmreport")

    # weitere globale Reports
    assert opts.glob_trusted_hitreport == str(
        rdir / "global_trusted_hits_summary.hmmreport"
    )
    assert opts.glob_intermediate_hitreport == str(
        rdir / "global_intermediate_hits_summary.hmmreport"
    )
    assert opts.csb_output_file == str(
        rdir / "Collinear_syntenic_blocks" / "Csb_output.txt"
    )


def test_prepare_result_space_preserves_existing_glob_report(tmp_path: Path):
    # Arrange: existierendes globales Report vorgeben
    preexisting = tmp_path / "pre.cat_hmmreport"
    preexisting.write_text("# dummy\n", encoding="utf-8")

    opts = Hmsss(
        location=str(tmp_path),
        result_files_directory=str(tmp_path / "results"),
        new_project=True,
        project_name="keepglob",
        glob_report=str(preexisting),  # bereits existierende Datei
    )

    # Act
    prepare_result_space(opts, project=opts.project_name)

    # Assert: glob_report wurde NICHT überschrieben
    assert opts.glob_report == str(preexisting)
    # Der Rest wurde dennoch korrekt aufgebaut
    rdir = Path(opts.result_files_directory)
    assert rdir.parent == tmp_path / "results"
    assert rdir.name.endswith("_keepglob")
    assert (rdir / "Filtered_hits").is_dir()


def test_prepare_result_space_custom_result_dir_when_not_new(tmp_path: Path):
    # Arrange: eigener Zielordner, new_project=False
    custom_root = tmp_path / "custom_results"
    opts = Hmsss(
        location=str(tmp_path),
        result_files_directory=str(custom_root),
        new_project=False,
        project_name="abc",
    )

    # Act
    prepare_result_space(opts, project=opts.project_name)

    # Assert: Projekt wurde unter custom_root angelegt
    rdir = Path(opts.result_files_directory)
    assert rdir.parent == custom_root
    assert rdir.name.endswith("_abc")

    # Felder/Verzeichnisse vorhanden
    assert opts.cross_check_directory == str(rdir / "Filtered_hits")
    assert (rdir / "Filtered_hits").is_dir()
    assert opts.glob_report == str(rdir / "global_report.cat_hmmreport")
