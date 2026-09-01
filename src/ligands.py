"""Converts ligand SMILES to 3D for binding"""


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import csv


METALS = {"Zn", "Ca", "Na", "K", "Mg", "Fe", "Cu", "Mn", "Al", "Ga",
          "As", "Se", "Xe", "Ag", "Au", "Pt", "Li", "Ba", "Bi", "Sb"
          "He", "Ne", "Ar", "Kr", "Xe", "Rn",}


def smiles_to_3d(smiles, out_path):
    """Build a 3D structure from a SMILES string and write as SDF"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    embed_molecule = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_molecule == -1:
        return None
    AllChem.MMFFOptimizeMolecule(mol)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()
    return out_path


def prepare_library(csv_path, out_dir):
    """Generate 3D structures for every SMILES in the CSV. Returns (n_ok, n_failed, n_skipped)."""
    batching_counter = 0
    failure_counter = 0
    skipped_counter = 0
    with open(csv_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if not is_dockable(row["smiles"]):
                skipped_counter += 1
                continue
            out_path = out_dir / f"{row['chembl_id']}.sdf"
            result = smiles_to_3d(row["smiles"], out_path)
            batching_counter += 1
            if batching_counter % 100 == 0:
                print(f"Batching count: {batching_counter}")
            if result is None:
                failure_counter += 1

    return (batching_counter - failure_counter, failure_counter, skipped_counter)



def is_dockable(smiles, max_mw=600.0):
    """True if the SMILES is a single, metal-free, reasonably sized organic molecule."""
    if "." in smiles:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALS:
            return False
    if Descriptors.MolWt(mol) > max_mw:
        return False
    return True