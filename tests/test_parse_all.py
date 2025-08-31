# tests/test_show_options.py
from __future__ import annotations

import os
import pathlib

from hmsss.cli.parse import _project_root_from_this_file
from hmsss.cli.paths import PROJECT_DIR


def test_project_root_from_this_file_returns_absolute_path():
    root = _project_root_from_this_file()
    print(root)
    # Rückgabewert ist ein String
    assert isinstance(root, str)
    # Rückgabewert ist ein absoluter Pfad
    assert os.path.isabs(root)
    assert root == PROJECT_DIR


def test_project_root_from_this_file_points_three_levels_above():
    here = pathlib.Path(__file__).resolve()
    # Der Pfad zur parse.py im Projekt
    parse_py = here.parents[2] / "src" / "hmsss" / "cli" / "parse.py"
    assert parse_py.exists(), "parse.py nicht gefunden, Testumgebung falsch"

    expected_root = parse_py.parent.parent.parent.resolve()
    returned_root = pathlib.Path(_project_root_from_this_file()).resolve()

    assert returned_root == expected_root


def test_project_root_from_this_file_exists_as_directory():
    root = pathlib.Path(_project_root_from_this_file())
    assert root.exists() and root.is_dir()
