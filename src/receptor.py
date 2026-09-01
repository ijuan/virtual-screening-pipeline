from pathlib import Path
import subprocess

def extract_chain(pdb_path, chain, out_path):
    """Writes only the ATOM records for `chain` from `pdb_path` to `out_path`.

    Drops HETATM records (water, ions, glycans, ligands) and all other chains.
    Returns the output path.
    """
    matching_list = []

    pdb_line = pdb_path.read_text().splitlines()

    for line in pdb_line:
        if (line.startswith("ATOM")) and (line[21] == chain):
            matching_list.append(line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out_file:
        out_file.write("\n".join(matching_list))
    return out_path


def extract_ligand(pdb_path, ligand_key, out_path):
    matching_list = []
    
    pdb_line = pdb_path.read_text().splitlines()

    for line in pdb_line:
        if (line.startswith("HETATM")) and ((line[17:20].strip(), line[21], 
                                           line[22:26].strip()) == ligand_key):
            matching_list.append(line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out_file:
        out_file.write("\n".join(matching_list))
    return out_path


def to_pdbqt(pdb_path, out_path, rigid=True):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["obabel", str(pdb_path), "-O", str(out_path), "-p", "7.4"]

    if rigid:
        cmd.append("-xr")

    subprocess.run(cmd, capture_output=True, check=True)
    return out_path


