"""Featurize all SMILES into Morgan's Fingerprints for XGBoost to run through"""


import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import csv
from rdkit.Chem import rdFingerprintGenerator

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def featurize(smiles, n_bits=2048, radius=2):
    """Turn a SMILES string into a numeric feature vector. Returns None on invalid input."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    arr = _MORGAN_GEN.GetFingerprintAsNumPy(mol)
    desc = np.array([
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumRotatableBonds(mol),
        mol.GetNumHeavyAtoms(),
    ])
    return np.concatenate([arr, desc])


def build_dataset(scores_csv, library_csv):
    """Featurize every scored molecule. Returns (X, y, ids)."""
    X = []
    y = []
    ids = []

    with open(library_csv) as f:
        reader = csv.DictReader(f)
        smiles_by_id = {row["chembl_id"]: row["smiles"] for row in reader}

    with open(scores_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row["score"]:
                continue
            featurization_smiles = featurize(smiles_by_id[row["chembl_id"]])
            if featurization_smiles is None:
                continue
            X.append(featurization_smiles)
            y.append(float(row["score"]))
            ids.append(row["chembl_id"])

    return np.array(X), np.array(y), ids


