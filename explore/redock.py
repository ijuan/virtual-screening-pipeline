from pathlib import Path

import src.structure as structure
import src.dock as dock
import src.validate as validate
from explore.prep_target import prepare_target


def redock(pdb_id, chain, receptor_dir, out_path, exhaustiveness=8):
    out_path = Path(out_path)
    key, box = structure.get_box_from_pdb(pdb_id)
    receptor_dir = Path(receptor_dir)
    receptor = receptor_dir / f"{pdb_id}_{chain}.pdbqt"
    ref_ligand = receptor_dir / f"{key[0]}_ref.pdbqt"

    center, size = box
    score = dock.dock_one(receptor, ref_ligand, center, size, out_path, exhaustiveness)

    docked_coords = validate.read_coords(out_path, first_model_only=True)
    ref_coords = validate.read_coords(ref_ligand)

    return score, validate.rmsd(docked_coords, ref_coords)


if __name__ == "__main__":
    print(redock("3ERT", "A", "data/receptor", "data/results/3ert/OHT_redock.pdbqt"))