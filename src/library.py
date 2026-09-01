import requests
import csv


def fetch_chembl_approved(out_path, limit=1000):
    """Fetch FDA-approved small molecules from ChEMBL, write to CSV."""
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4"
    offset = 0
    rows = []

    while True:
        params = {
            "max_phase": 4,
            "molecule_type": "Small molecule",
            "limit": limit,
            "offset": offset,
        }
        library_csv = requests.get(url, params=params)
        if library_csv.status_code != 200:
            raise ValueError("Could not fetch FDA library.")
        
        fda_data = library_csv.json()

        for molecule in fda_data["molecules"]:
            structures = molecule["molecule_structures"]
            if structures is None:
                continue
            canonical_smiles = structures["canonical_smiles"]
            chembl_id = molecule["molecule_chembl_id"]
            rows.append((chembl_id, canonical_smiles))
        if fda_data["page_meta"]["next"] is None:
            break
        offset += limit

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["chembl_id", "smiles"])
        writer.writerows(rows)

    return out_path