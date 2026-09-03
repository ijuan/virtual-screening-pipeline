"""Setup docking for all ligands"""

import csv
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def dock_one(receptor_path, ligand_path, center, size, out_path, exhaustiveness=8):
    """Dock one ligand into the receptor. Returns the best affinity in kcal/mol."""

    cmd = [
        "vina",
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", str(out_path),
        "--exhaustiveness", str(exhaustiveness),
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(cmd, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        return None

    for line in out_path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def dock_library(receptor_path, ligand_dir, center, size, out_dir, results_csv, workers=None):
    """Dock every ligand in ligand_dir. Writes (chembl_id, score) rows to CSV."""
    ligands = sorted(ligand_dir.glob("*.pdbqt"))
    results = []

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {}
        for lig in ligands:
            out_path = out_dir / lig.name
            fut = pool.submit(dock_one, receptor_path, lig, center, size, out_path)
            futures[fut] = lig.stem

        progress_counter = 0
        for future in as_completed(futures):
            chembl_id = futures[future]
            score = future.result()
            results.append((chembl_id, score))
            progress_counter += 1
            if progress_counter % 50 == 0:
                print(progress_counter)

    with open(results_csv, "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["chembl_id", "score"])
        csv_writer.writerows(results)
    return results_csv
    