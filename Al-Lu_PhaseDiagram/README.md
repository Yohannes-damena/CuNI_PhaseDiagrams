# Al–Lu Phase Diagram (PyCalphad)

This repository contains a Jupyter notebook that uses PyCalphad to:
- Load an Al–Lu `.tdb` database
- Generate the Al–Lu binary phase diagram
- Identify and annotate invariant reactions (eutectic / peritectic candidates)
- Compute equilibrium phase fractions for Al-30Lu vs temperature (°C)
- Provide interpretation guidance for alloy design

## 1) Environment setup

You need Python 3.10+.

Install dependencies (recommended in a virtual environment):

```bash
pip install -r requirements.txt
```

If you use Conda:

```bash
conda create -n al_lu python=3.11 -y
conda activate al_lu
pip install -r requirements.txt
```

## 2) Files
- `Al-LU.tdb`: Thermodynamic database (already included here)
- `Al_Lu_PhaseDiagram.ipynb`: Jupyter notebook (created in this repo)

## 3) Run the notebook

```bash
jupyter notebook
```
Then open `Al_Lu_PhaseDiagram.ipynb` and run all cells top-to-bottom.

## 4) Notes
- The notebook will attempt to consider all condensed phases present in the database for accurate topology. If performance is slow, you can reduce the phase list to `['LIQUID', 'FCC_A1', 'HCP_A3']` within the notebook.
- Invariant reaction classification is heuristic (based on just-above/just-below temperatures at the detected point). It works well in practice, but if any points are ambiguous, confirm by inspecting nearby phase fields in the diagram.


