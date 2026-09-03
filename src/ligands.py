"""Converts ligand SMILES to 3D for binding"""


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import csv
from src.receptor import to_pdbqt
import subprocess


ALLOWED_ELEMENTS = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"}


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
    """True if the SMILES is a single, made only of common elements, reasonably sized organic molecule."""
    if "." in smiles:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ALLOWED_ELEMENTS:
            return False
    if Descriptors.MolWt(mol) > max_mw:
        return False
    return True


def convert_library(sdf_dir, out_dir):
    """Convert every SDF in sdf_dir to a flexible PDBQT. Returns (n_ok, n_failed)."""
    failed = 0
    ok = 0
    for sdf_path in sdf_dir.glob("*.sdf"):
        try:
            out_path = out_dir / f"{sdf_path.stem}.pdbqt"
            to_pdbqt(sdf_path, out_path, rigid=False)
            ok += 1
        except subprocess.CalledProcessError:
            failed += 1
        if (ok + failed) % 100 == 0:
            print(f"Converted: {ok + failed}")
    return (ok, failed)