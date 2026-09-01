"""Find RMSD box is from crytal"""

def read_coords(path, first_model_only=False):
    path_lines = path.read_text().splitlines()
    coords_list = []
    for line in path_lines:
        if first_model_only and line.startswith("ENDMDL"):
            break
        if line[77:79].strip() in ("H", "HD"):
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords_list.append((x, y, z))
    return coords_list