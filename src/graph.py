from typing import TypedDict
import src.structure as structure
import explore.prep_target as target
import explore.redock as redock
import explore.screen as screen
import src.featurization as featurization
import numpy as np
import explore.active_learning as active_learning
from pathlib import Path


class PipelineState(TypedDict):
    pdb_id: str
    chain: str
    pdb_path: str
    ligand_key: tuple
    box: list
    receptor_pdbqt: str
    ref_ligand_pdbqt: str
    redock_score: float
    rmsd: float
    scores_csv: str
    features_path: str
    history: list
    receptor_dir: str      # "data/receptor"
    results_dir: str       # "data/results/3ert"


def find_box(state):
    key, box = structure.get_box_from_pdb(state["pdb_id"])
    return {"ligand_key": key, "box": box}


def prep_receptor(state):
    receptor_pdbqt, ref_ligand_pdbqt = target.prepare_target(
        state["pdb_id"], state["chain"], state["ligand_key"], state["receptor_dir"]
    )
    return {"receptor_pdbqt": receptor_pdbqt, "ref_ligand_pdbqt": ref_ligand_pdbqt}


def redock_node(state):
    out_path = Path(state["results_dir"]) / "redock.pdbqt"
    score, rmsd = redock.redock(
        state["receptor_pdbqt"], state["ref_ligand_pdbqt"], state["box"], out_path
    )
    return {"redock_score": score, "rmsd": rmsd}


def screen_node(state):
    out_csv = Path(state["results_dir"]) / "scores.csv"
    pose_dir = Path(state["results_dir"]) / "poses"
    screen.screen(
        state["receptor_pdbqt"], state["ligand_dir"], state["box"],
        pose_dir, out_csv, workers=state.get("workers", 8),
    )
    return {"scores_csv": str(out_csv)}


def featurize_node(state):
    X, y, ids = featurization.build_dataset(state["scores_csv"], state["library_csv"])
    features_path = Path(state["results_dir"]) / "features.npz"
    np.savez(features_path, X=X, y=y, ids=np.array(ids))
    return {"features_path": str(features_path)}


def active_learning_node(state):
    data = np.load(state["features_path"], allow_pickle=True)
    history = active_learning.run_from_arrays(
        data["X"], data["y"], list(data["ids"]),
        active_learning.select_uncertainty,
    )
    return {"history": history}

