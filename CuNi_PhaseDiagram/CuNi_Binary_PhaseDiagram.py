# --- Import Setup ---
from pycalphad import Database, variables as v, binplot, equilibrium
import matplotlib.pyplot as plt
import numpy as np


# --- Load the Database ---

db = Database("data/CuNi_db.tdb")

# ---Define the components and phases ---
# --- VA for vaccum ---
components = ['CU', 'NI', 'VA']
phases = ['LIQUID', 'FCC_A1']

# --- Binary Phase Diagram ---

fig, ax = plt.subplots(figsize=(8,6))
binplot(db, components, phases, x=v.x('NI'), y=v.T,
        plot_kwargs={'color': 'royalblue', 'lw':2})
ax.set_xlabel('Mole Fraction Ni')
ax.set_ylabel('Temprature(K)')
ax.set_title('Cu-Ni Binary Phase Diagram')
ax.grid(True, alpha=0.3)
plt.show()

# -- Equilibrum phase Fractions for Cu-40Ni

T_range = np.linspace(400, 1800, 200)
conds = {v.X('NI'): 0.40, v.T: T_range, v.P: 101325}
eq = equilibrium(db, components, phases, conds)

fig, ax = plt.subplots(figsize=(8,6))
for phase in phases:
     ax.plot(eq.T - 273.15, eq.NP.where(eq.phase == phase).sum(dim='vertex'),
             label=phase)

ax.set_xlabel('Temperature (°C)')
ax.set_ylabel('Phase Fraction')
ax.set_title('Phase Fractions vs Temperature for Cu-40Ni')
ax.legend()
plt.show()


