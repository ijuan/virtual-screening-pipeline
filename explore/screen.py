from pathlib import Path

import src.structure as structure
import src.dock as dock


def screen(receptor_pdbqt, ligand_dir, box, pose_dir, out_csv, workers=8):
    center, size = box
    return dock.dock_library(
        Path(receptor_pdbqt), Path(ligand_dir), center, size,
        Path(pose_dir), Path(out_csv), workers=workers,
    )


if __name__ == "__main__":
    print(screen(
        "3ERT", "A",
        "data/receptor",
        "data/ligands/pdbqt",
        "data/results/3ert/poses",
        "data/results/3ert/scores.csv",
    ))