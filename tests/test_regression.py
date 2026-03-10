import pathlib

from helena import run


def test_AR2plus(snaptolshot):
    import os

    os.chdir("tests/")

    output_file = pathlib.Path("data/TECPlot2D/2DPlots_Data/AR2+.csv")

    output_file.unlink(missing_ok=True)

    run(argv=[])

    assert output_file.exists()

    with output_file.open("r") as f:
        output_lines = f.readlines()

    assert len(output_lines) > 0

    assert snaptolshot == output_lines

    pathlib.Path("data/meshnodes.dat").unlink()
    pathlib.Path("data/TECPlot2D/2DPlots_Data/AR2+.csv").unlink()
    pathlib.Path("data/TECPlot2D/2DPlot_AR2+.png").unlink()
    pathlib.Path("data/TECPlot2D/2DPlots_Data").rmdir()
    pathlib.Path("data/TECPlot2D").rmdir()

    os.chdir("../")
