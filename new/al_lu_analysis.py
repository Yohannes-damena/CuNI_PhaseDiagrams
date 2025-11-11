import sys
import math
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.plot.binary import binplot
import xarray as xr


def load_system(db_path: str) -> Tuple[Database, List[str], List[str], Dict]:
	"""
	Load the Al-Lu database and define system settings with restricted phases.
	"""
	db = Database(db_path)
	components = ['AL', 'LU', 'VA']
	phases = ['LIQUID', 'FCC_A1', 'HCP_A3']
	conditions = {v.P: 101325}  # Pa
	return db, components, phases, conditions


def make_binplot(
	db: Database,
	components: List[str],
	phases: List[str],
	t_min: float = 500.0,
	t_max: float = 3000.0,
	n_t: int = 121,
	save_path: str = "al_lu_binary.png"
) -> Tuple[plt.Figure, plt.Axes]:
	"""
	Generate the Al-Lu binary phase diagram with restricted phases.
	"""
	fig = plt.figure(figsize=(8, 6), dpi=150)
	# binplot will calculate the phase diagram projection using the provided phases
	# Some pycalphad versions expect the composition variable to be the first X() in conditions.
	# Use X(AL) grid and relabel/invert axis to display X(Lu).
	conditions = {
		v.P: 101325,
		v.T: np.linspace(t_min, t_max, n_t),
		v.X('AL'): np.linspace(0.0, 1.0, 201)
	}
	binplot(db, components, phases, conditions)
	ax = plt.gca()
	# Invert x-axis so that it increases with X(Lu) = 1 - X(AL)
	ax.set_xlim(1.0, 0.0)
	ax.set_xlabel('X(Lu)')
	ax.set_ylabel('Temperature (K)')
	ax.set_title('Al–Lu Binary Phase Diagram (Restricted: LIQUID, FCC_A1, HCP_A3)')
	ax.grid(True, alpha=0.2, linestyle='--')
	fig.tight_layout()
	fig.savefig(save_path)
	return fig, ax


def compute_eq_grid(
	db: Database,
	components: List[str],
	phases: List[str],
	conditions_base: Dict,
	x_points: np.ndarray,
	t_points: np.ndarray
) -> xr.Dataset:
	"""
	Compute an equilibrium grid over composition and temperature.
	Returns an xarray Dataset with phase amounts under eq['NP'].
	"""
	conds = dict(conditions_base)
	conds.update({
		v.X('LU'): x_points,
		v.T: t_points
	})
	ds = equilibrium(db, components, phases, conds, output='NP')
	return ds


def _phase_fractions_from_eq(eq: xr.Dataset, tol: float = 1e-8) -> Dict[str, xr.DataArray]:
	"""
	Extract per-phase fractions from eq['NP'] by normalizing with total moles.
	Returns a dict mapping phase name to DataArray over the remaining condition dims.
	"""
	if 'NP' not in eq.variables or 'Phase' not in eq.variables:
		raise ValueError("Equilibrium dataset missing required variables 'NP' or 'Phase'.")
	np_da = eq['NP']  # dims typically include ('N','P','T','X_*','vertex')
	phase_da = eq['Phase']  # same dims as NP, each 'vertex' labeled with a phase name
	# Total moles (sum over vertices)
	if 'vertex' in np_da.dims:
		total_moles = np_da.sum(dim='vertex')
	else:
		total_moles = np_da
	total_moles_safe = xr.where(total_moles > tol, total_moles, np.nan)
	# Determine unique phase names present in dataset
	unique_phases = set([str(x) for x in np.unique(phase_da.values)])
	phase_fractions: Dict[str, xr.DataArray] = {}
	for ph in unique_phases:
		# Indicator for this phase in each vertex
		indicator = xr.where(phase_da.astype(str) == ph, 1.0, 0.0)
		amount_this = (np_da * indicator)
		if 'vertex' in amount_this.dims:
			amount_this = amount_this.sum(dim='vertex')
		frac = (amount_this / total_moles_safe).fillna(0.0)
		phase_fractions[ph] = frac
	return phase_fractions


def detect_invariant_clusters(
	phase_fractions: Dict[str, xr.DataArray],
	tol: float = 1e-3
) -> List[Tuple[int, int]]:
	"""
	Detect grid points where LIQUID, FCC_A1, and HCP_A3 are all present (approximate three-phase).
	Returns list of (i_x, i_t) indices representing cluster centers.
	Clusters neighboring True cells into a single point per cluster via BFS.
	"""
	required = ['LIQUID', 'FCC_A1', 'HCP_A3']
	for r in required:
		if r not in phase_fractions:
			return []
	mask = (phase_fractions['LIQUID'] > tol) & (phase_fractions['FCC_A1'] > tol) & (phase_fractions['HCP_A3'] > tol)
	# Determine index names for X and T (handle different dimension orders)
	dims = list(mask.dims)
	# Find closest matches for 'X' and 'T' in dims
	x_dim = next((d for d in dims if d.upper().startswith('X')), None)
	t_dim = next((d for d in dims if d.upper().startswith('T')), None)
	if x_dim is None or t_dim is None:
		return []
	mask_np = mask.transpose(x_dim, t_dim).values
	nx, nt = mask_np.shape
	visited = np.zeros_like(mask_np, dtype=bool)
	cluster_centers: List[Tuple[int, int]] = []
	neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
	for i in range(nx):
		for j in range(nt):
			if not mask_np[i, j] or visited[i, j]:
				continue
			# BFS to get a cluster
			queue = [(i, j)]
			visited[i, j] = True
			points = [(i, j)]
			while queue:
				ci, cj = queue.pop(0)
				for di, dj in neighbors:
					ni, nj = ci + di, cj + dj
					if 0 <= ni < nx and 0 <= nj < nt and not visited[ni, nj] and mask_np[ni, nj]:
						visited[ni, nj] = True
						queue.append((ni, nj))
						points.append((ni, nj))
			# Center of mass (average indices)
			pi = int(round(np.mean([p[0] for p in points])))
			pj = int(round(np.mean([p[1] for p in points])))
			cluster_centers.append((pi, pj))
	return cluster_centers


def classify_invariant(
	db: Database,
	components: List[str],
	phases: List[str],
	conditions_base: Dict,
	x_val: float,
	t_val: float,
	dt: float = 5.0,
	tol: float = 1e-3
) -> str:
	"""
	Classify an invariant around (x_val, t_val) as 'Eutectic' or 'Peritectic' based on local equilibria.
	"""
	def stable_set_at(T_probe: float) -> Set[str]:
		ds = equilibrium(db, components, phases, {
			**conditions_base,
			v.X('LU'): x_val,
			v.T: T_probe
		}, output='NP')
		pf = _phase_fractions_from_eq(ds, tol=tol)
		actives = {ph for ph, arr in pf.items() if float(arr) > tol}
		return actives

	above = stable_set_at(t_val + dt)
	below = stable_set_at(t_val - dt)
	# Eutectic: above is primarily Liquid (optionally tiny solid), below has two solids and no liquid
	if ('LIQUID' in above) and (('FCC_A1' not in above) or ('HCP_A3' not in above)) and \
		('LIQUID' not in below) and (('FCC_A1' in below) and ('HCP_A3' in below)):
		return 'Eutectic'
	# Peritectic: above has Liquid + one solid; below has the other solid only (or dominant)
	if ('LIQUID' in above) and (('FCC_A1' in above) ^ ('HCP_A3' in above)) and \
		('LIQUID' not in below) and (('FCC_A1' in below) ^ ('HCP_A3' in below)):
		return 'Peritectic'
	# Fallback classification
	return 'Invariant'


def annotate_invariants_and_solids(
	ax: plt.Axes,
	db: Database,
	components: List[str],
	phases: List[str],
	conditions_base: Dict,
	x_points: np.ndarray,
	t_points: np.ndarray,
	eq_grid: xr.Dataset,
	marker_size: int = 60
) -> List[Tuple[float, float, str]]:
	"""
	Annotate detected invariants (eutectic/peritectic) and mark representative solid-state transitions.
	Returns list of (x, T, label) for invariants.
	"""
	phase_fracs = _phase_fractions_from_eq(eq_grid)
	centers = detect_invariant_clusters(phase_fracs)
	invariant_markers: List[Tuple[float, float, str]] = []
	if len(centers) > 0:
		# Determine dim order for mapping indices back to values
		sample_arr = next(iter(phase_fracs.values()))
		dims = list(sample_arr.dims)
		x_dim = next((d for d in dims if d.upper().startswith('X')), None)
		t_dim = next((d for d in dims if d.upper().startswith('T')), None)
		for (ix, it) in centers:
			x_val = float(x_points[ix])
			t_val = float(t_points[it])
			label = classify_invariant(db, components, phases, conditions_base, x_val, t_val)
			ax.scatter(x_val, t_val, s=marker_size, marker='*', c='red', edgecolor='k', zorder=5)
			ax.text(x_val, t_val + 30, label, ha='center', va='bottom', fontsize=8, color='red')
			invariant_markers.append((x_val, t_val, label))

	# Solid-state transitions: mark where dominant solid switches between FCC_A1 and HCP_A3 in regions with negligible liquid
	liquid = phase_fracs.get('LIQUID')
	fcc = phase_fracs.get('FCC_A1')
	hcp = phase_fracs.get('HCP_A3')
	if liquid is not None and fcc is not None and hcp is not None:
		liquid_np = liquid.transpose(...).values
		fcc_np = fcc.transpose(...).values
		hcp_np = hcp.transpose(...).values
		# Ensure arrays are X x T
		# Try to interpret dims
		dims = list(liquid.dims)
		x_dim = next((d for d in dims if d.upper().startswith('X')), None)
		t_dim = next((d for d in dims if d.upper().startswith('T')), None)
		if x_dim is not None and t_dim is not None:
			liquid_np = liquid.transpose(x_dim, t_dim, ...).values
			fcc_np = fcc.transpose(x_dim, t_dim, ...).values
			hcp_np = hcp.transpose(x_dim, t_dim, ...).values
			nx, nt = liquid_np.shape
			for i in range(nx):
				# Consider temperatures where liquid is negligible
				for j in range(1, nt):
					if liquid_np[i, j] < 1e-4 and liquid_np[i, j - 1] < 1e-4:
						dominant_prev = 'FCC_A1' if fcc_np[i, j - 1] >= hcp_np[i, j - 1] else 'HCP_A3'
						dominant_curr = 'FCC_A1' if fcc_np[i, j] >= hcp_np[i, j] else 'HCP_A3'
						if dominant_prev != dominant_curr:
							ax.scatter(float(x_points[i]), float(t_points[j]), s=20, marker='s', c='black', zorder=4)
							# Only mark first transition per composition to avoid clutter
							break
	return invariant_markers


def run_30lu_fraction_plot(
	db: Database,
	components: List[str],
	phases: List[str],
	conditions_base: Dict,
	t_min: float = 400.0,
	t_max: float = 2000.0,
	n_t: int = 801,
	x_lu: float = 0.30,
	save_path: str = "al_lu_phase_fractions_0p30.png"
) -> Dict[str, np.ndarray]:
	"""
	Compute and plot phase fractions versus temperature at X(LU)=x_lu.
	Returns a mapping of phase -> (T[K], fraction) arrays for downstream analysis.
	"""
	t_points = np.linspace(t_min, t_max, n_t)
	ds = equilibrium(db, components, phases, {
		**conditions_base,
		v.X('LU'): x_lu,
		v.T: t_points
	}, output='NP')
	phase_fracs = _phase_fractions_from_eq(ds)
	# Prepare plot
	t_k = t_points
	t_c = t_k - 273.15
	plt.figure(figsize=(8, 5), dpi=150)
	for ph in ['LIQUID', 'FCC_A1', 'HCP_A3']:
		if ph in phase_fracs:
			# Squeeze to 1D over T
			arr = phase_fracs[ph]
			# Drop all dims except T
			if 'T' in arr.dims:
				y = np.asarray(arr.transpose('T', ...).values).ravel()
			else:
				# Try to find dim like 'T'
				t_dim = next((d for d in arr.dims if d.upper().startswith('T')), None)
				y = np.asarray(arr.transpose(t_dim, ...).values).ravel() if t_dim else np.zeros_like(t_c)
			plt.plot(t_c, y, label=ph)
	# Styling
	plt.xlabel('Temperature (°C)')
	plt.ylabel('Phase fraction')
	plt.title(f'Al–Lu Phase Fractions at X(Lu) = {x_lu:.2f} (Restricted phases)')
	plt.xlim(t_c.min(), t_c.max())
	plt.ylim(0, 1)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.legend()
	plt.tight_layout()
	plt.savefig(save_path)
	plt.close()
	# Return arrays for analysis
	result = {}
	for ph in ['LIQUID', 'FCC_A1', 'HCP_A3']:
		if ph in phase_fracs:
			arr = phase_fracs[ph]
			if 'T' in arr.dims:
				y = np.asarray(arr.transpose('T', ...).values).ravel()
			else:
				t_dim = next((d for d in arr.dims if d.upper().startswith('T')), None)
				y = np.asarray(arr.transpose(t_dim, ...).values).ravel() if t_dim else np.zeros_like(t_k)
		else:
			y = np.zeros_like(t_k)
		result[ph] = np.vstack([t_k, y])
	return result


def identify_transformations_at_x(result: Dict[str, np.ndarray], tol: float = 1e-3) -> Dict[str, List[Tuple[float, float]]]:
	"""
	From phase fractions vs T, identify approximate transformation temperatures.
	Returns dict with keys 'eutectic', 'peritectic' mapping to list of (T[K], T[°C]) tuples.
	"""
	t_k = result['LIQUID'][0]
	liq = result['LIQUID'][1]
	fcc = result['FCC_A1'][1] if 'FCC_A1' in result else np.zeros_like(t_k)
	hcp = result['HCP_A3'][1] if 'HCP_A3' in result else np.zeros_like(t_k)
	out = {'eutectic': [], 'peritectic': []}
	for i in range(1, len(t_k) - 1):
		above = i + 1
		below = i - 1
		active_above = set()
		active_below = set()
		if liq[above] > tol:
			active_above.add('LIQUID')
		if fcc[above] > tol:
			active_above.add('FCC_A1')
		if hcp[above] > tol:
			active_above.add('HCP_A3')
		if liq[below] > tol:
			active_below.add('LIQUID')
		if fcc[below] > tol:
			active_below.add('FCC_A1')
		if hcp[below] > tol:
			active_below.add('HCP_A3')
		# Eutectic signature
		if ('LIQUID' in active_above) and (len(active_above - {'LIQUID'}) <= 0) and \
			('LIQUID' not in active_below) and (('FCC_A1' in active_below) and ('HCP_A3' in active_below)):
			out['eutectic'].append((float(t_k[i]), float(t_k[i] - 273.15)))
		# Peritectic signature
		if ('LIQUID' in active_above) and (('FCC_A1' in active_above) ^ ('HCP_A3' in active_above)) and \
			('LIQUID' not in active_below) and (('FCC_A1' in active_below) ^ ('HCP_A3' in active_below)):
			out['peritectic'].append((float(t_k[i]), float(t_k[i] - 273.15)))
	return out


def print_analysis(x_lu: float, invariants: List[Tuple[float, float, str]], transforms_30lu: Dict[str, List[Tuple[float, float]]]) -> None:
	"""
	Print a concise analysis and interpretation based on results.
	"""
	print("=== Analysis and Interpretation ===")
	if invariants:
		print("Detected invariant points (approximate):")
		for x_val, t_val, label in invariants:
			print(f" - {label} near X(Lu)≈{x_val:.3f}, T≈{t_val:.1f} K ({t_val - 273.15:.1f} °C)")
	else:
		print("No clear three-phase invariants detected within the restricted phase set.")
	print()
	print(f"Solidification path for Al-{int(x_lu*100)}Lu (X(Lu)={x_lu:.2f}):")
	is_eut = len(transforms_30lu.get('eutectic', [])) > 0
	is_peri = len(transforms_30lu.get('peritectic', [])) > 0
	if is_eut:
		tk, tc = transforms_30lu['eutectic'][0]
		print(f" - Eutectic-like event at ≈ {tk:.1f} K ({tc:.1f} °C).")
	if is_peri:
		tk, tc = transforms_30lu['peritectic'][0]
		print(f" - Peritectic-like event at ≈ {tk:.1f} K ({tc:.1f} °C).")
	if not (is_eut or is_peri):
		print(" - No strong eutectic/peritectic signature within restricted phases; likely single-solid formation followed by solid-state transition.")
	print(" - Phases on cooling appear/disappear according to the fraction curves; observe where LIQUID vanishes and solid phases emerge.")
	print()
	print("Comparison to Cu–Ni (isomorphous):")
	print(" - Cu–Ni exhibits a simple FCC isomorphous system with no eutectic/peritectic.")
	print(" - In contrast, Al–Lu (restricted here to FCC/HCP/LIQUID) can show two-solid fields and potential invariant reactions.")
	print()
	print("Role of rare-earth (Lu) in stability/kinetics:")
	print(" - Rare-earth elements often stabilize intermetallic/close-packed structures, shifting transition temperatures.")
	print(" - They can increase liquidus/solidus separation and influence nucleation barriers, affecting transformation kinetics.")
	print(" - Practical implication: control cooling rates to tailor microstructures (e.g., avoid coarse eutectics or promote fine lamellae).")
	print()
	print("Alloy design/processing implications:")
	print(" - Knowledge of invariants guides casting/heat-treatment windows and segregation control.")
	print(" - Solid-state FCC↔HCP transition markers suggest potential for phase-transformation strengthening or texture control.")
	print(" - For precise paths, include all database phases in a follow-up study; here we restricted to LIQUID/FCC_A1/HCP_A3 as requested.")


def main() -> None:
	# 1) Load database and system definition
	db, components, phases, base_conds = load_system('Al-Lu.tdb')
	# 2) Binary phase diagram
	x_points = np.linspace(0.0, 1.0, 121)
	t_points = np.linspace(500.0, 3000.0, 121)
	try:
		fig, ax = make_binplot(db, components, phases, t_min=500.0, t_max=3000.0, n_t=161, save_path='al_lu_binary.png')
		# Equilibrium grid for annotations
		eq_grid = compute_eq_grid(db, components, phases, base_conds, x_points, t_points)
	except Exception as e:
		# Fallback: build a phase-region map from equilibrium
		print(f"binplot failed, using fallback map: {e}")
		eq_grid = compute_eq_grid(db, components, phases, base_conds, x_points, t_points)
		phase_fracs = _phase_fractions_from_eq(eq_grid)
		liq = phase_fracs.get('LIQUID')
		fcc = phase_fracs.get('FCC_A1')
		hcp = phase_fracs.get('HCP_A3')
		# Build Z map: 0=LIQUID, 1=FCC, 2=HCP
		if (liq is None) or (fcc is None) or (hcp is None):
			raise RuntimeError("Missing phase fractions for fallback plotting.")
		# Ensure arrays are X x T
		# Drop singleton N/P dims if present
		for dim_to_drop in ['N', 'P']:
			if dim_to_drop in liq.dims:
				liq = liq.isel({dim_to_drop: 0})
			if dim_to_drop in fcc.dims:
				fcc = fcc.isel({dim_to_drop: 0})
			if dim_to_drop in hcp.dims:
				hcp = hcp.isel({dim_to_drop: 0})
		dims = list(liq.dims)
		x_dim = next((d for d in dims if d.upper().startswith('X')), None)
		t_dim = next((d for d in dims if d.upper().startswith('T')), None)
		liq_np = liq.transpose(x_dim, t_dim, ...).values
		fcc_np = fcc.transpose(x_dim, t_dim, ...).values
		hcp_np = hcp.transpose(x_dim, t_dim, ...).values
		nx, nt = liq_np.shape
		Z = np.zeros_like(liq_np, dtype=int)
		for i in range(nx):
			for j in range(nt):
				if liq_np[i, j] >= 1e-3:
					Z[i, j] = 0
				else:
					Z[i, j] = 1 if fcc_np[i, j] >= hcp_np[i, j] else 2
		Xg, Tg = np.meshgrid(x_points, t_points, indexing='ij')
		fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
		cmap = plt.get_cmap('viridis', 3)
		im = ax.pcolormesh(Xg, Tg, Z, cmap=cmap, shading='nearest')
		cbar = plt.colorbar(im, ax=ax, ticks=[0.33, 1, 1.66], fraction=0.046, pad=0.04)
		cbar.ax.set_yticklabels(['LIQUID', 'FCC_A1', 'HCP_A3'])
		ax.set_xlabel('X(Lu)')
		ax.set_ylabel('Temperature (K)')
		ax.set_title('Al–Lu Binary (Fallback map; restricted phases)')
		ax.grid(True, alpha=0.2, linestyle='--')
		ax.set_xlim(0.0, 1.0)
		ax.set_ylim(t_points.min(), t_points.max())
		fig.tight_layout()
		fig.savefig('al_lu_binary.png')
	invariants = annotate_invariants_and_solids(ax, db, components, phases, base_conds, x_points, t_points, eq_grid)
	fig.savefig('al_lu_binary.png')  # Save again with annotations
	plt.close(fig)
	# 3) Phase fractions at Al-30Lu, plot in °C
	result = run_30lu_fraction_plot(db, components, phases, base_conds, t_min=400.0, t_max=2000.0, n_t=801, x_lu=0.30, save_path='al_lu_phase_fractions_0p30.png')
	# 4) Identify transformations for X(Lu)=0.30 and print analysis
	transforms_30lu = identify_transformations_at_x(result)
	print_analysis(0.30, invariants, transforms_30lu)


if __name__ == '__main__':
	main()


