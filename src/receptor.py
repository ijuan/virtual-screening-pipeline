from pathlib import Path


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

    with open(out_path, "w") as out_file:
        out_file.write("\n".join(matching_list))
    return matching_list