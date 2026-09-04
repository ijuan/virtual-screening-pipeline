"""Sort Molecules into similar groups. Assign groups to train or test 80/20"""


from rdkit.Chem.Scaffolds import MurckoScaffold
import collections


def get_scaffold(smiles):
    """Return the Bemis-Murcko scaffold as a SMILES string."""
    try:
       return MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles)
    except Exception:
        return None


def scaffold_split(ids, smiles_by_id, test_frac=0.2):
    """Split IDs into (train_ids, test_ids) so no scaffold appears in both."""
    train_ids = []
    test_ids = []
    scaffold_dict = collections.defaultdict(list)
    for chembl_id in ids:
        scaffold = get_scaffold(smiles_by_id[chembl_id])
        if scaffold is None:
            continue
        scaffold_dict[scaffold].append(chembl_id)

    
    sorted_scaffold = sorted(scaffold_dict.values(), key=len, reverse=True)
    for group in sorted_scaffold:
        if len(test_ids) < test_frac * len(ids):
            test_ids.extend(group)
        else:
            train_ids.extend(group)
    return (train_ids, test_ids)
