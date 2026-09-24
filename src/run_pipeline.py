from src.graph import build_graph

if __name__ == "__main__":
    app = build_graph()
    result = app.invoke({
        "pdb_id": "3ERT",
        "chain": "A",
        "receptor_dir": "data/receptor",
        "results_dir": "data/results/3ert",
        "ligand_dir": "data/ligands/pdbqt",
        "library_csv": "data/library/fda.csv",
        "workers": 8,
    })
    print(result["rmsd"], result.get("history"))