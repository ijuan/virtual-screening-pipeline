from pathlib import Path

import src.structure as structure
import src.dock as dock


def screen(pdb_id, chain, receptor_dir, ligand_dir, out_csv,
           workers=8, exhaustiveness=8, padding=3.0):
    receptor_dir = Path(receptor_dir)
    out_csv = Path(out_csv)

    key, box = structure.get_box_from_pdb(pdb_id)
    center, size = box
    receptor = receptor_dir / f"{pdb_id}_{chain}.pdbqt"

    return dock.dock_library(
        receptor, Path(ligand_dir), center, size, out_csv,
        workers=workers, exhaustiveness=exhaustiveness,
    )


if __name__ == "__main__":
    print(screen(
        "3ERT", "A",
        "data/receptor",
        "data/ligands/pdbqt",
        "data/results/3ert/scores.csv",
    ))