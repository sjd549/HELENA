import os
import shutil
import pathlib
import pytest


@pytest.fixture(scope="session")
def input_directory():
    return pathlib.Path(__file__).parent / "input_files"


@pytest.fixture(scope="session")
def data_directory():
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_env(tmp_path_factory, data_directory):
    """
    Copy the expensive/shared input data once per test session.
    """
    root = tmp_path_factory.mktemp("helena_env")
    shutil.copytree(data_directory, root / "data")
    return root


@pytest.fixture
def work_env(tmp_path, test_env):
    """
    Give each test its own isolated copy of the seeded environment.
    """
    workdir = tmp_path / "work"
    shutil.copytree(test_env / "data", workdir / "data")
    return workdir


@pytest.fixture
def run_in_tmp_dir(work_env, monkeypatch):
    """
    Temporarily change the directory to the working directory.
    """
    monkeypatch.chdir(work_env)
