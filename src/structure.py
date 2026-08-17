import pathlib
import requests

def fetch_pdb(pdb_id, data_dir="data/structures"):
    pdb_id = pdb_id.upper()
    target_path = pathlib.Path(data_dir) / f"{pdb_id}.pdb"
    if target_path.exists():
        return target_path
    target_path.parent.mkdir(parents=True,exist_ok=True)
    
    response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")

    # Checks for 404 error before writing into file
    if response.status_code != 200:
        raise ValueError(f"Could not fetch {pdb_id}")

    with open(target_path, "w") as f:
        f.write(response.text)

    return target_path

