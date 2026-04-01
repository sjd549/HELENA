import pytest
from pydantic import ValidationError

from helena import run


@pytest.mark.parametrize(
    "input_file",
    ["input_tecplot2D.toml"],
)
def test_AR2plus(input_file, input_directory, snaptolshot, work_env_with_chdir):
    run(argv=[str(input_directory / input_file)])

    output_file = work_env_with_chdir / "data/TECPlot2D/2DPlots_Data/AR2+.csv"
    assert output_file.exists()

    output_lines = output_file.read_text().splitlines(keepends=True)
    assert output_lines
    assert snaptolshot == output_lines


@pytest.mark.usefixtures("work_env_with_chdir")
def test_input(input_directory):
    with pytest.raises(ValidationError):
        run(argv=[str(input_directory / "input_wrong_type.toml")])
