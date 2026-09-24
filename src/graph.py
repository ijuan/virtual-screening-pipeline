from typing import TypedDict
import src.structure as structure
import explore.prep_target as target
import explore.redock as redock
import explore.screen as screen
import src.featurization as featurization
import numpy as np
import explore.active_learning as active_learning
from pathlib import Path
from langgraph.graph import StateGraph, END


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
    status: str
    ligand_dir: str
    library_csv: str
    workers: int



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


RMSD_THRESHOLD = 2.0


def gate_on_rmsd(state):
    """Conditional edge: only screen if the redock validated."""
    return "screen" if state["rmsd"] <= RMSD_THRESHOLD else "failed"


def failed_node(state):
    return {"status": f"redock failed validation: {state['rmsd']:.2f} A RMSD"}


def build_graph():
    g = StateGraph(PipelineState)

    g.add_node("find_box", find_box)
    g.add_node("prep_receptor", prep_receptor)
    g.add_node("redock", redock_node)
    g.add_node("screen", screen_node)
    g.add_node("featurize", featurize_node)
    g.add_node("active_learning", active_learning_node)
    g.add_node("failed", failed_node)

    g.set_entry_point("find_box")
    g.add_edge("find_box", "prep_receptor")
    g.add_edge("prep_receptor", "redock")
    g.add_conditional_edges("redock", gate_on_rmsd)
    g.add_edge("screen", "featurize")
    g.add_edge("featurize", "active_learning")
    g.add_edge("active_learning", END)
    g.add_edge("failed", END)

    return g.compile()

