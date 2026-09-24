from pathlib import Path

import src.structure as structure
import src.dock as dock
import src.validate as validate
from explore.prep_target import prepare_target


def redock(receptor_pdbqt, ref_ligand_pdbqt, box, out_path, exhaustiveness=8):
    out_path = Path(out_path)
    center, size = box
    score = dock.dock_one(receptor_pdbqt, ref_ligand_pdbqt, center, size, out_path, exhaustiveness)
    docked = validate.read_coords(out_path, first_model_only=True)
    ref = validate.read_coords(ref_ligand_pdbqt)
    return score, validate.rmsd(docked, ref)


if __name__ == "__main__":
    print(redock("3ERT", "A", "data/receptor", "data/results/3ert/OHT_redock.pdbqt"))