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
            chain        = line[21]
            res_num      = line[22:26].strip()
            key = (residue_name, chain, res_num)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            parsed_hetatm.setdefault(key, []).append((x, y, z))
    return parsed_hetatm

# Set of residue names that aren't in the ligand
IGNORE = {"HOH", "GOL", "SO4", "PO4", "EDO", "PEG", "ACT",
          "ZN", "NA", "MG", "CL", "CA", "K", "MN", "FE",
          "NAG", "BMA", "MAN", "FUC", "GAL", "SIA", "GLC"}

# Returns which residues are the actual drug
def pick_ligand(hetatm: dict, min_atoms=10, ligand=None) -> str:
    if ligand is not None:
        key_matches = [k for k in hetatm if k[0] == ligand]
        if not key_matches:
            raise ValueError(f"{ligand} not found. Present: {sorted({k[0] for k in hetatm})}")
        return sorted(key_matches)[0]
    filtered_hetatm = []
    for key, coords in hetatm.items():
        residue_name , chain, res_num = key
        if len(coords) < min_atoms:
            continue
        if residue_name in IGNORE:
            continue
        filtered_hetatm.append(key)

    if not filtered_hetatm:
        raise ValueError(f"No ligand found - structure may be apo." 
                            f" HETATM residues present: {sorted({k[0] for k in hetatm})}")
        
    best_candidate = None
    for name in filtered_hetatm:
        if best_candidate is None or len(hetatm[name]) > len(hetatm[best_candidate]):
            best_candidate = name
    return best_candidate


def compute_box(coords:list, padding=8.0) -> list:
    x_coords = [c[0] for c in coords]
    y_coords = [c[1] for c in coords]
    z_coords = [c[2] for c in coords]

    x_avg = sum(x_coords) / len(x_coords)
    y_avg = sum(y_coords) / len(y_coords)
    z_avg = sum(z_coords) / len(z_coords)
    center = (x_avg, y_avg, z_avg)

    x_size = (max(x_coords) - min(x_coords)) + 2 * padding
    y_size = (max(y_coords) - min(y_coords)) + 2 * padding
    z_size = (max(z_coords) - min(z_coords)) + 2 * padding
    size = ( x_size, y_size, z_size)

    return [center, size]






if __name__ == "__main__":
   x = compute_box([(0.639, 22.832, 7.09), (1.119, 22.308, 8.239), (0.248, 22.065, 9.241), (2.415, 22.022, 8.41), (3.34, 22.233, 7.433)])
   print(x)