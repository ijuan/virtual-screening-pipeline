"""Normalizes inflated scores based off number of heavy atoms"""
from rdkit import Chem
import csv


def heavy_atom_count(smiles):
    """Number of non-hydrogen atoms in the molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol.GetNumHeavyAtoms()


def add_efficiency(scores_csv, library_csv, out_csv, min_atoms=15):
    """Join scores with SMILES, compute ligand efficiency, write sorted CSV."""
    list_heavy_atom_score = []
    with open(library_csv) as f:
        reader = csv.DictReader(f)
        smiles_by_id = {row["chembl_id"]: row["smiles"] for row in reader}
    with open(scores_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            chembl_id = row["chembl_id"]
            if not row["score"]:
                continue
            score = float(row["score"])
            n_atoms = heavy_atom_count(smiles_by_id[chembl_id])
            if n_atoms < min_atoms:
                continue
            efficiency = score / n_atoms
            list_heavy_atom_score.append((chembl_id, score, n_atoms, efficiency))

    list_heavy_atom_score.sort(key=lambda r: r[3])
    with open(out_csv, "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["chembl_id", "score", "n_atoms", "efficiency"])
        csv_writer.writerows(list_heavy_atom_score)
    return out_csv
