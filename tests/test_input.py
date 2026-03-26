from helena import run


def test_AR2plus(snaptolshot, input_directory, work_env_with_chdir):
    input_file = input_directory / "input1.toml"

    run(argv=[str(input_file)])

    output_file = work_env_with_chdir / "data/TECPlot2D/2DPlots_Data/AR2+.csv"
    assert output_file.exists()

    output_lines = output_file.read_text().splitlines(keepends=True)
    assert output_lines
    assert snaptolshot == output_lines
