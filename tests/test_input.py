import pathlib
import pytest
import shutil
import os

from helena import run


@pytest.fixture
def directory():
    return pathlib.Path(__file__).parent / "input_files"


def test_AR2plus(snaptolshot, directory):
    input_file = directory / "input1.toml"

    os.chdir("tests/")

    output_file = pathlib.Path("data/TECPlot2D/2DPlots_Data/AR2+.csv")

    output_file.unlink(missing_ok=True)

    run(argv=[str(input_file)])

    assert output_file.exists()

    with output_file.open("r") as f:
        output_lines = f.readlines()

    assert len(output_lines) > 0

    assert snaptolshot == output_lines

    pathlib.Path("data/meshnodes.dat").unlink()
    shutil.rmtree(pathlib.Path("data/TECPlot2D/"))

    os.chdir("../")
