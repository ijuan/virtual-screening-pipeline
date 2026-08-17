# Structure-Based Virtual Screening Pipeline

Given a protein structure (PDB ID), screen compound libraries and return a ranked
shortlist of candidate binders.

**Status:** in development. Environment setup complete; pipeline not yet implemented.

---

## Approach

Molecular docking (AutoDock Vina) scores how well a compound fits into a protein's
binding pocket. Docking is accurate but slow — roughly 10-60 seconds per compound on
a laptop, which makes exhaustive screening of a 400k-compound library infeasible.

This pipeline uses docking scores as training labels for a machine learning model,
which then predicts scores for undocked compounds. An active-learning loop uses those
predictions to decide which compounds are worth docking next, recovering most of the
top hits while docking only a fraction of the library.

### Planned funnel

```
PDB ID
  ↓  fetch structure, locate binding pocket from co-crystallized ligand
docking box
  ↓  FDA-approved library (~2,500) — docked exhaustively
  ↓  natural product library (COCONUT) — active-learning loop
docked results
  ↓  binding affinity threshold
  ↓  Lipinski rule-of-5 filter
ranked shortlist
```

**Binding affinity** is reported in kcal/mol. More negative = stronger binding.

### Scope decisions

- **Binding site from co-crystallized ligand.** The pocket center is derived from the
  coordinates of a ligand already bound in the PDB structure. This is reliable and
  deterministic, but requires a holo structure — apo structures (protein alone) are
  rejected rather than guessed at.
- **Vina retained over a pure-ML approach.** Docking provides physics-based scoring and,
  more importantly, generates training labels without wet-lab work.

---

## Setup

### 1. Create the environment

Requires conda (Miniforge or Anaconda).

```bash
conda env create -f environment.yml
conda activate docking
```

### 2. Verify

```bash
python -c "from rdkit import Chem; from vina import Vina; from meeko import MoleculePreparation; import pymol; print('all ok')"
```

If this prints `all ok`, the environment is working.

---

## Why conda instead of pip

Most of this stack is compiled C++ with Python bindings — RDKit, Open Babel, Vina,
PyMOL. These install cleanly from conda-forge and are unreliable via pip, which has
to fall back to building from source when no wheel matches the platform and Python
version.

**One exception: `meeko` is installed via pip**, declared in the `pip:` section of
`environment.yml`. Meeko is pure Python, so conda offers no advantage, and the
conda-forge build is years out of date — it imports `rdkit.six`, a module removed
from RDKit long ago, and fails on import. PyPI has the current release.

Because pip does not resolve conda-managed dependencies, meeko's own requirements
(`scipy`, `gemmi`) are listed explicitly as conda dependencies rather than left to pip.

**Python is pinned to 3.12.** Newer releases are ahead of the meeko and Vina binding
release cycles.

---

## Layout

```
src/          pipeline modules
explore/      scratch work and evidence — not imported by src/
data/         libraries, structures, results (gitignored)
tests/
main.py       entry point
```

`data/` is gitignored, as are structure files (`.pdb`, `.pdbqt`), results (`.csv`),
and PyMOL sessions (`.pse`). Docking generates these in bulk and they do not belong
in version control.

---

## Tools

| Tool | Role |
|---|---|
| RDKit | SMILES parsing, 3D conformer generation, descriptors, Lipinski |
| Meeko | RDKit molecule → PDBQT (atom types, charges, rotatable bonds) |
| AutoDock Vina | docking and scoring |
| Biopython | fetching and parsing PDB structures |
| PyMOL | visual verification of box placement and poses |
| scikit-learn / XGBoost | affinity prediction model, active-learning loop |