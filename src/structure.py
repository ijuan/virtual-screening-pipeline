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


# Pull out HETATM lines
def parse_hetatms(pdb_path) -> dict:
    parsed_hetatm = {}
    pdb_line = pdb_path.read_text().splitlines()

    for line in pdb_line:
        if line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            parsed_hetatm.setdefault(residue_name,[]).append((x, y, z))
    return parsed_hetatm

# Set of residue names that aren't in the ligand
IGNORE = {"HOH", "GOL", "SO4", "PO4", "EDO", "PEG", "ACT",
          "ZN", "NA", "MG", "CL", "CA", "K", "MN", "FE"}

# Returns which residues are the actual drug
def pick_ligand(hetatm: dict, min_atoms=10) -> str:
    filtered_hetatm = []
    for name, coords in hetatm.items():
        if len(coords) < min_atoms:
            continue
        if name in IGNORE:
            continue
        filtered_hetatm.append(name)

    if not filtered_hetatm:
        raise ValueError(f"No ligand found - structure may be apo." 
                            f" HETATM residues present: {list(hetatm)}")
        
    best = None
    for name in filtered_hetatm:
        if best is None or len(hetatm[name]) > len(hetatm[best]):
            best = name
    return best

        

    







