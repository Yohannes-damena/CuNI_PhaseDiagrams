# Al–Lu Phase Diagram (PyCalphad)

This repository contains a CLI Python script that uses PyCalphad to:
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
- `al_lu_analysis.py`: CLI script (binary diagram + Al-30Lu equilibrium)

## 3) Run the script

```bash
# Default usage (expects Al-LU.tdb in the same folder)
python al_lu_analysis.py

# Specify options
python al_lu_analysis.py --db Al-LU.tdb --xlu 0.30 --tmin 400 --tmax 2000 --binplot-tmin 500 --binplot-tmax 3000

# Restrict phases to a minimal core set (faster)
python al_lu_analysis.py --phases core
```

Outputs are saved under `outputs/`:
- `al_lu_binary.png`: Binary phase diagram with heuristic E/P annotations
- `al30lu_phasefractions.png`: Al-30Lu phase fractions vs Temperature (°C)
- `al30lu_transformations.txt`: Detected invariant temperatures (if any)

## 4) Notes
- The script will consider all phases present in the database by default. If performance is slow, pass `--phases core` to limit to `['LIQUID', 'FCC_A1', 'HCP_A3']`.
- Invariant reaction classification is heuristic (based on just-above/just-below temperatures at the detected point). It works well in practice, but if any points are ambiguous, confirm by inspecting nearby phase fields in the diagram.
- If your `Al-LU.tdb` is a unary (pure-elements) database (e.g., SGTE Unary), it lacks Al–Lu interaction parameters and intermetallics; real eutectic/peritectic features will not be reproduced. Use an assessed Al–Lu database for realistic results.


