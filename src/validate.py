"""Find RMSD box is from crytal"""

import math


def read_coords(path, first_model_only=False):
    path_lines = path.read_text().splitlines()
    coords_dict = {}
    for line in path_lines:
        if first_model_only and line.startswith("ENDMDL"):
            break
        if line[77:79].strip() in ("H", "HD"):
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords_dict[line[12:16].strip()] = (x, y, z)
    return coords_dict


def rmsd(coords_a, coords_b):
    """Calculate rmsd from ligand to actual"""
    squared_list = []
    shared = set(coords_a) & set(coords_b)
    if len(coords_a) != len(coords_b):
        raise ValueError
    for name in shared:
        xa, ya, za = coords_a[name]
        xb, yb, zb = coords_b[name]
        squared_dist = (xa-xb)**2 + (ya-yb)**2 + (za-zb)**2
        squared_list.append(squared_dist)
    mean = sum(squared_list) / len(squared_list)
    return math.sqrt(mean) # RMSD