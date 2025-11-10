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

fig, ax = plt.subplots(figsize=(10,8))
binplot(db, components, phases, {v.X('NI'): (0, 1, 0.01), v.T: (400, 1800, 10), v.P: 101325},
        x=v.X('NI'), y=v.T, ax=ax,
        plot_kwargs={'color': 'royalblue', 'lw':3})

# Enhanced axis labels and formatting
ax.set_xlabel('Mole Fraction of Ni (x$_{Ni}$)', fontsize=14, fontweight='bold')
ax.set_ylabel('Temperature (K)', fontsize=14, fontweight='bold')
ax.set_title('Cu-Ni Binary Phase Diagram', fontsize=16, fontweight='bold', pad=20)

# Add axis labels for pure components
ax.text(0.02, 400, 'Cu', fontsize=12, fontweight='bold', ha='center')
ax.text(0.98, 400, 'Ni', fontsize=12, fontweight='bold', ha='center')

# Add phase region labels
ax.text(0.2, 1200, 'Liquid (L)', fontsize=12, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
ax.text(0.7, 600, 'FCC-A1\n(Solid Solution)', fontsize=12, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7))
ax.text(0.5, 1000, 'L + FCC-A1\n(Two-Phase Region)', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))

# Add annotations for key temperatures
ax.axhline(y=1358, color='red', linestyle='--', alpha=0.7, linewidth=2)
ax.text(0.05, 1380, 'T$_{m,Cu}$ = 1358 K', fontsize=10, color='red', fontweight='bold')

ax.axhline(y=1728, color='red', linestyle='--', alpha=0.7, linewidth=2)
ax.text(0.05, 1750, 'T$_{m,Ni}$ = 1728 K', fontsize=10, color='red', fontweight='bold')

# Add grid and formatting
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_xlim(0, 1)
ax.set_ylim(400, 1800)

# Add minor ticks
ax.minorticks_on()
ax.tick_params(which='minor', length=3)
ax.tick_params(which='major', length=6, width=1.5)

plt.tight_layout()
plt.savefig('CuNi_Binary_PhaseDiagram.png', dpi=400, bbox_inches='tight')
plt.close()

# =============================================================================
# EQUILIBRIUM PHASE FRACTION SIMULATION
# =============================================================================
# Fixed composition: Cu-40Ni (X(Ni) = 0.40)
# Temperature range: 400 K to 1800 K
# Convert temperature to Celsius for clarity

print("="*60)
print("EQUILIBRIUM PHASE FRACTION SIMULATION")
print("="*60)
print(f"Composition: Cu-40Ni (X(Ni) = 0.40)")
print(f"Temperature range: 400 K to 1800 K")
print("="*60)

# Define temperature range with higher resolution for better accuracy
T_range = np.linspace(400, 1800, 300)  # Increased resolution
conds = {v.X('NI'): 0.40, v.T: T_range, v.P: 101325}

print("Computing equilibrium phase fractions...")
eq = equilibrium(db, components, phases, conds)
print("Calculation completed!")

# Create the phase fraction plot
fig, ax = plt.subplots(figsize=(12, 8))

# Plot phase fractions for each phase
phase_colors = {'LIQUID': 'red', 'FCC_A1': 'blue'}
phase_labels = {'LIQUID': 'Liquid (L)', 'FCC_A1': 'FCC-A1 (Solid Solution)'}

for phase in phases:
    phase_fraction = eq.NP.where(eq.Phase == phase).sum(dim='vertex')
    ax.plot(eq.T - 273.15, phase_fraction.squeeze(),
            label=phase_labels[phase],
            linewidth=4,
            marker='o',
            markersize=3,
            color=phase_colors[phase],
            alpha=0.8)

# Enhanced formatting and annotations
ax.set_xlabel('Temperature (°C)', fontsize=16, fontweight='bold')
ax.set_ylabel('Phase Fraction', fontsize=16, fontweight='bold')
ax.set_title('Equilibrium Phase Fractions vs Temperature\nCu-40Ni Alloy (X(Ni) = 0.40)',
             fontsize=18, fontweight='bold', pad=25)

# Add composition annotation
ax.text(0.02, 0.95, f'Composition: Cu-40Ni\nX(Ni) = 0.40',
        transform=ax.transAxes, fontsize=12, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))

# Add horizontal reference lines
ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.7, linewidth=2)
ax.text(ax.get_xlim()[1]*0.7, 0.52, '50% Phase Fraction',
        fontsize=11, color='gray', fontweight='bold')

ax.axhline(y=1.0, color='black', linestyle='-', alpha=0.3, linewidth=1)
ax.axhline(y=0.0, color='black', linestyle='-', alpha=0.3, linewidth=1)

# Calculate and annotate key transition temperatures
liquid_phase = eq.NP.where(eq.Phase == 'LIQUID').sum(dim='vertex').squeeze()
fcc_phase = eq.NP.where(eq.Phase == 'FCC_A1').sum(dim='vertex').squeeze()

# Find solidus temperature (where liquid first appears)
solidus_temp_idx = np.where(liquid_phase > 0.01)[0]  # 1% threshold
if len(solidus_temp_idx) > 0:
    solidus_temp = eq.T[solidus_temp_idx[0]] - 273.15
    ax.axvline(x=solidus_temp, color='red', linestyle='--', alpha=0.8, linewidth=3)
    ax.text(solidus_temp + 30, 0.85, f'SOLIDUS\n{solidus_temp:.0f}°C',
            fontsize=12, color='red', fontweight='bold', ha='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    print(f"Solidus temperature: {solidus_temp:.1f}°C")

# Find liquidus temperature (where solid disappears)
liquidus_temp_idx = np.where(fcc_phase < 0.99)[0]  # 99% threshold
if len(liquidus_temp_idx) > 0:
    liquidus_temp = eq.T[liquidus_temp_idx[-1]] - 273.15
    ax.axvline(x=liquidus_temp, color='blue', linestyle='--', alpha=0.8, linewidth=3)
    ax.text(liquidus_temp - 30, 0.85, f'LIQUIDUS\n{liquidus_temp:.0f}°C',
            fontsize=12, color='blue', fontweight='bold', ha='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    print(f"Liquidus temperature: {liquidus_temp:.1f}°C")

# Calculate melting range
if len(solidus_temp_idx) > 0 and len(liquidus_temp_idx) > 0:
    melting_range = liquidus_temp - solidus_temp
    ax.text(0.5, 0.15, f'Melting Range: {melting_range:.0f}°C',
            transform=ax.transAxes, fontsize=12, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.8))
    print(f"Melting range: {melting_range:.1f}°C")

# Enhanced legend and formatting
ax.legend(fontsize=14, loc='center right', framealpha=0.9,
          fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_ylim(-0.05, 1.05)
ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1])

# Add minor ticks and formatting
ax.minorticks_on()
ax.tick_params(which='minor', length=4)
ax.tick_params(which='major', length=8, width=2, labelsize=12)

# Add phase region annotations
ax.text(0.15, 0.5, 'Two-Phase\nRegion', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7),
        ha='center', va='center')
ax.text(0.85, 0.5, 'Single-Phase\nRegions', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7),
        ha='center', va='center')

plt.tight_layout()
plt.savefig('CuNi_PhaseFractions.png', dpi=400, bbox_inches='tight')
plt.close()

print("="*60)
print("SIMULATION COMPLETED")
print("="*60)


