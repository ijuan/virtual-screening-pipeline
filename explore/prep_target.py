from pathlib import Path

import src.structure as structure
import src.receptor as receptor


def prepare_target(pdb_id, chain, out_dir):
    out_dir = Path(out_dir)
    pdb = structure.fetch_pdb(pdb_id)
    key, box = structure.get_box_from_pdb(pdb_id)

    chain_pdb = out_dir / f"{pdb_id}_{chain}.pdb"
    receptor.extract_chain(pdb, chain, chain_pdb)
    out_pdbqt = out_dir / f"{pdb_id}_{chain}.pdbqt"
    receptor.to_pdbqt(chain_pdb, out_pdbqt)

    ligand_pdb = out_dir / f"{key[0]}_ref.pdb"
    receptor.extract_ligand(pdb, key, ligand_pdb)
    ligand_pdbqt = out_dir / f"{key[0]}_ref.pdbqt"
    receptor.to_pdbqt(ligand_pdb, ligand_pdbqt, rigid=False)

    return key, box


if __name__ == "__main__":
    print(prepare_target("3ERT", "A", "data/receptor"))