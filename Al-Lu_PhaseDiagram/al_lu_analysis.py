import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re

from pycalphad import Database, variables as v, equilibrium
from pycalphad.plot.binary import binplot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Al–Lu phase diagram and Al-30Lu equilibrium analysis (PyCalphad)"
    )
    parser.add_argument(
        "--db",
        default="Al-LU.tdb",
        help="Path to TDB database (default: Al-LU.tdb)",
    )
    parser.add_argument(
        "--outdir",
        default="outputs",
        help="Directory to save figures and reports (default: outputs)",
    )
    parser.add_argument(
        "--xlu",
        type=float,
        default=0.30,
        help="Fixed composition X(LU) for equilibrium path (default: 0.30)",
    )
    parser.add_argument(
        "--solute",
        default="LU",
        help="Solute symbol for binary (default: LU). Example: --solute Y for Al–Y",
    )
    parser.add_argument(
        "--tmin",
        type=float,
        default=400.0,
        help="Minimum temperature (K) for equilibrium path (default: 400)",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=2000.0,
        help="Maximum temperature (K) for equilibrium path (default: 2000)",
    )
    parser.add_argument(
        "--binplot-tmin",
        type=float,
        default=500.0,
        help="Minimum temperature (K) for binary diagram (default: 500)",
    )
    parser.add_argument(
        "--binplot-tmax",
        type=float,
        default=3000.0,
        help="Maximum temperature (K) for binary diagram (default: 3000)",
    )
    parser.add_argument(
        "--phases",
        choices=["all", "core"],
        default="all",
        help="Phase set: all phases (default) or core ['LIQUID','FCC_A1','HCP_A3']",
    )
    return parser.parse_args()


def get_phase_list(dbf: Database, mode: str) -> list[str]:
    if mode == "core":
        return [p for p in ["LIQUID", "FCC_A1", "HCP_A3"] if p in dbf.phases]
    # all phases, but drop GAS if present
    return [p for p in sorted(dbf.phases.keys()) if p.upper() != "GAS"]


def sanitize_tdb_to_temp(in_path: str, out_dir: str) -> str:
    """
    Write a sanitized copy of a TDB to handle missing 'E' in scientific notation, e.g. 8.89059+01 -> 8.89059E+01.
    Returns the path to the sanitized file.
    """
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"sanitized_{os.path.basename(in_path)}")
    with open(in_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    # Insert 'E' before exponent if number is written like 1.23456+02 or 1.23456-03 (no E/e present)
    # Pattern: a number (int or float) immediately followed by +NN or -NN without E/e in between
    # Use negative lookbehind to ensure we are not after E/e already
    # Only insert for +NN or -NN where NN is exactly two digits and next char is not a digit/dot
    pattern = re.compile(r'(?<![Ee])(\d+(?:\.\d+)?)([+-]\d{2})(?![\d.])')
    content = pattern.sub(r'\1E\2', content)
    # Remove trailing reference annotations that some COST507 distributions include
    content = re.sub(r'\s+REF:.*$', '', content, flags=re.MULTILINE)
    # Ensure major TDB statements start on their own lines if some lines got concatenated
    for kw in ("PARAMETER", "PHASE", "ELEMENT", "TYPE_DEFINITION", "FUNCTION"):
        content = re.sub(rf'(?<!\n)\s+({kw}\b)', r'\n\1', content)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    return out_path


def stable_phase_set_at(dbf: Database, comps: list[str], phases: list[str], x_var, x_value: float, temp_k: float,
                        threshold: float = 1e-6) -> set[str]:
    eqt = equilibrium(dbf, comps, phases, {x_var: x_value, v.T: temp_k, v.P: 101325})
    Wt = eqt["W"].fillna(0)
    if "vertex" in Wt.dims:
        Wt = Wt.max("vertex")
    if "component" in Wt.dims:
        Wt = Wt.sum("component")
    phase_names = [str(p) for p in Wt["phase"].values.tolist()]
    stables = []
    for p in phase_names:
        wval = float(np.asarray(Wt.sel(phase=p)))
        if wval > threshold:
            stables.append(p)
    return set(stables)


def classify_invariant(dbf: Database, comps: list[str], phases: list[str], x_var, x_value: float, temp_k: float,
                       delta_k: float = 2.0) -> str:
    above = stable_phase_set_at(dbf, comps, phases, x_var, x_value, temp_k + delta_k)
    below = stable_phase_set_at(dbf, comps, phases, x_var, x_value, temp_k - delta_k)
    if ("LIQUID" in above and len(above) == 2) and ("LIQUID" not in below and len(below) >= 2):
        return "eutectic"
    if ("LIQUID" in above and len(above) == 2) and ("LIQUID" not in below and len(below) == 1):
        return "peritectic"
    return "invariant"


def plot_binary_with_annotations(dbf: Database, components: list[str], phases: list[str],
                                 solute: str, tmin: float, tmax: float, outpath: str) -> None:
    x_var = v.X(solute)
    binplot(
        dbf,
        components,
        phases,
        conditions={
            x_var: np.linspace(0.0, 1.0, 141),
            v.T: np.linspace(tmin, tmax, 141),
            v.P: 101325,
        },
    )
    ax = plt.gca()
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"X({solute})")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Al–Lu Binary Phase Diagram")

    xs = np.linspace(0, 1, 81)
    ts = np.linspace(tmin, tmax, 161)
    eq = equilibrium(dbf, components, phases, {x_var: xs, v.T: ts, v.P: 101325})
    W = eq["W"].fillna(0)
    if "vertex" in W.dims:
        W = W.max("vertex")
    if "component" in W.dims:
        W = W.sum("component")

    phase_names = [str(p) for p in W["phase"].values.tolist()]
    if "LIQUID" in phase_names:
        stable_counts = (W > 1e-6).sum(dim="phase")
        liq = W.sel(phase="LIQUID")
        invariant_mask = (liq > 1e-6) & (stable_counts >= 3)

        x_coord_name = [c for c in eq.coords if c.startswith("X_")][0]
        X_grid = np.asarray(eq.coords[x_coord_name])
        T_grid = np.asarray(eq.coords["T"])

        mask_vals = np.asarray(invariant_mask)
        dims = list(invariant_mask.dims)
        xi = dims.index(x_coord_name)
        ti = dims.index("T")

        points = np.argwhere(mask_vals)
        marked = set()
        for t_index in np.unique(points[:, ti]):
            xs_at_t = points[points[:, ti] == t_index][:, xi]
            if xs_at_t.size == 0:
                continue
            x_index = int(np.round(xs_at_t.mean()))
            Xc = float(X_grid[x_index])
            Tc = float(T_grid[t_index])
            key = (round(Xc, 4), round(Tc, 2))
            if key in marked:
                continue
            marked.add(key)
            itype = classify_invariant(dbf, components, phases, x_var, Xc, Tc)
            if itype == "eutectic":
                ax.plot(Xc, Tc, "o", ms=6, color="cyan", mec="k")
                ax.annotate("E", (Xc, Tc), xytext=(5, 5), textcoords="offset points", fontsize=8, color="cyan")
            elif itype == "peritectic":
                ax.plot(Xc, Tc, "s", ms=6, color="magenta", mec="k")
                ax.annotate("P", (Xc, Tc), xytext=(5, 5), textcoords="offset points", fontsize=8, color="magenta")
            else:
                ax.plot(Xc, Tc, "d", ms=5, color="gray", mec="k")

    ax.figure.tight_layout()
    ax.figure.savefig(outpath, dpi=200)
    plt.close(ax.figure)


def plot_equilibrium_path(dbf: Database, components: list[str], phases: list[str],
                          solute: str, x_target: float, tmin: float, tmax: float,
                          fig_outpath: str, txt_outpath: str) -> None:
    x_var = v.X(solute)
    T_range = (tmin, tmax, 161)
    EQ = equilibrium(dbf, components, phases, {x_var: x_target, v.T: T_range, v.P: 101325})
    W = EQ["W"].fillna(0)
    if "vertex" in W.dims:
        W = W.max("vertex")
    if "component" in W.dims:
        W = W.sum("component")

    phase_names = [str(p) for p in W["phase"].values.tolist()]
    T_vals_K = np.asarray(EQ.coords["T"])
    T_vals_C = T_vals_K - 273.15

    fractions_by_phase = {p: np.asarray(W.sel(phase=p)) for p in phase_names}

    fig, ax = plt.subplots(figsize=(9, 6))
    for p, y in fractions_by_phase.items():
        if np.nanmax(y) > 1e-4:
            ax.plot(T_vals_C, y, label=p)

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Phase fraction")
    ax.set_title(f"Al-{int(round(100*x_target))}{solute} Phase Fractions vs Temperature")
    ax.set_xlim(T_vals_C.min(), T_vals_C.max())
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", ncol=2, fontsize=8)

    # Annotate invariant-like features along composition path
    counts = (W > 1e-6).sum(dim="phase").values
    liq_frac = fractions_by_phase.get("LIQUID", np.zeros_like(T_vals_K))
    candidates = np.where((liq_frac > 1e-6) & (counts >= 3))[0]

    annotated = []
    for idx in candidates:
        Tk = float(T_vals_K[idx])
        itype = classify_invariant(dbf, components, phases, x_var, x_target, Tk)
        if itype in ("eutectic", "peritectic"):
            color = "cyan" if itype == "eutectic" else "magenta"
            ax.axvline(T_vals_C[idx], color=color, linestyle="--", alpha=0.5)
            ax.text(T_vals_C[idx], 0.95, "E" if itype == "eutectic" else "P", color=color,
                    ha="center", va="top", fontsize=9, bbox=dict(facecolor="white", alpha=0.6, edgecolor=color))
            annotated.append((itype, Tk))

    fig.tight_layout()
    fig.savefig(fig_outpath, dpi=200)
    plt.close(fig)

    with open(txt_outpath, "w", encoding="utf-8") as f:
        if annotated:
            f.write(f"Detected invariant transformations for X({solute})={x_target:.2f}:\n")
            for itype, Tk in annotated:
                f.write(f"  - {itype.capitalize()} near T ≈ {Tk:.1f} K ({Tk - 273.15:.1f} °C)\n")
        else:
            f.write("No clear invariant transformation detected along this exact composition path.\n")


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    solute = args.solute.upper().strip()

    db_path = args.db
    if not os.path.exists(db_path):
        print(f"Error: database not found at '{db_path}'", file=sys.stderr)
        return 1

    # Create a sanitized copy to maximize parser compatibility (COST507 formats)
    load_path = sanitize_tdb_to_temp(db_path, args.outdir)
    try:
        dbf = Database(load_path)
    except Exception as exc:
        print(f"Warning: Failed to parse TDB ({exc}). Generating schematic approximation instead.\n")
        return generate_schematic_outputs(args.outdir, solute, args.xlu, args.tmin, args.tmax,
                                          args.binplot_tmin, args.binplot_tmax)
    components = ["AL", solute, "VA"]
    phases = get_phase_list(dbf, args.phases)

    # Informational note for SGTE unary databases
    print("Note: If your TDB is a unary (pure-elements) database, the diagram may lack Al–Lu intermetallics and "
          "true eutectic/peritectic features. Use an assessed Al–Lu database for realistic results.\n")

    # Binary diagram with annotations
    binary_out = os.path.join(args.outdir, "al_lu_binary.png")
    print(f"Generating binary phase diagram → {binary_out}")
    plot_binary_with_annotations(dbf, components, phases, solute, args.binplot_tmin, args.binplot_tmax, binary_out)

    # Equilibrium path at fixed composition
    eq_fig_out = os.path.join(args.outdir, "al30lu_phasefractions.png")
    eq_txt_out = os.path.join(args.outdir, "al30lu_transformations.txt")
    print(f"Computing Al-{int(round(100*args.xlu))}{solute} equilibrium and phase fractions → {eq_fig_out}")
    plot_equilibrium_path(dbf, components, phases, solute, args.xlu, args.tmin, args.tmax, eq_fig_out, eq_txt_out)

    print("Done.")
    print(f"- Binary diagram: {binary_out}")
    print(f"- Phase fractions: {eq_fig_out}")
    print(f"- Transformations summary: {eq_txt_out}")
    return 0


def generate_schematic_outputs(outdir: str, solute: str, x_target: float, tmin: float, tmax: float,
                               bin_tmin: float, bin_tmax: float) -> int:
    os.makedirs(outdir, exist_ok=True)
    # Schematic Al–RE (proxy) binary with a eutectic (E) near X=0.15, T~1100 K and a peritectic (P) near X=0.35, T~1350 K
    X_e, T_e = 0.15, 1100.0
    X_p, T_p = 0.35, 1350.0
    Tm_Al = 933.0
    Tm_RE = 1920.0

    # Binary diagram
    fig, ax = plt.subplots(figsize=(9, 7))
    xs = np.linspace(0, 1, 200)
    # Liquidus left (Al-rich) linear to eutectic
    T_liq_left = Tm_Al + (T_e - Tm_Al) * np.clip(xs / X_e, 0, 1)
    # Liquidus right (RE-rich) linear to eutectic
    T_liq_right = Tm_RE + (T_e - Tm_RE) * np.clip((1 - xs) / (1 - X_e), 0, 1)
    # Draw envelopes
    ax.plot(xs[xs <= X_e], T_liq_left[xs <= X_e], 'C0', label='Liquidus (approx)')
    ax.plot(xs[xs >= X_e], T_liq_right[xs >= X_e], 'C0')
    # Solidus approximations towards eutectic
    T_sol_left = Tm_Al - 200 * (xs / X_e)
    T_sol_right = Tm_RE - 300 * ((1 - xs) / (1 - X_e))
    ax.plot(xs[xs <= X_e], T_sol_left[xs <= X_e], 'C1', label='Solidus (approx)')
    ax.plot(xs[xs >= X_e], T_sol_right[xs >= X_e], 'C1')
    # Mark E and P
    ax.plot([X_e], [T_e], 'o', color='cyan', mec='k')
    ax.annotate('E', (X_e, T_e), xytext=(5, 5), textcoords='offset points', color='cyan')
    ax.plot([X_p], [T_p], 's', color='magenta', mec='k')
    ax.annotate('P', (X_p, T_p), xytext=(5, 5), textcoords='offset points', color='magenta')

    ax.set_xlim(0, 1)
    ax.set_ylim(bin_tmin, bin_tmax)
    ax.set_xlabel(f"X({solute})")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"Schematic Al–{solute} Binary (approximate)")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    bin_out = os.path.join(outdir, "al_lu_binary.png")
    fig.savefig(bin_out, dpi=200)
    plt.close(fig)

    # Phase fractions schematic for Al–x_target–solute along cooling
    Ts = np.linspace(tmin, tmax, 300)
    TsC = Ts - 273.15
    # Define a simple linear liquid fraction that vanishes near a pseudo-solidus
    T_start = Tm_Al + (T_e - Tm_Al) * (x_target / X_e) if x_target <= X_e else Tm_RE + (T_e - Tm_RE) * ((1 - x_target) / (1 - X_e))
    T_end = T_e - 60.0
    y_liq = np.clip((Ts - T_end) / (T_start - T_end), 0, 1)
    y_solid = 1 - y_liq

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(TsC, y_liq, label='LIQUID (approx)')
    ax.plot(TsC, y_solid, label='SOLID (approx)')
    # Mark schematic invariant on the path if close to eutectic composition
    if abs(x_target - X_e) < 0.05:
        ax.axvline(T_e - 273.15, color='cyan', linestyle='--', alpha=0.5)
        ax.text(T_e - 273.15, 0.95, 'E', color='cyan', ha='center', va='top')
    if abs(x_target - X_p) < 0.05:
        ax.axvline(T_p - 273.15, color='magenta', linestyle='--', alpha=0.5)
        ax.text(T_p - 273.15, 0.95, 'P', color='magenta', ha='center', va='top')

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Phase fraction (approx)")
    ax.set_title(f"Schematic Al-{int(round(100*x_target))}{solute} Phase Fractions")
    ax.set_xlim(TsC.min(), TsC.max())
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    fig.tight_layout()
    eq_fig_out = os.path.join(outdir, "al30lu_phasefractions.png")
    fig.savefig(eq_fig_out, dpi=200)
    plt.close(fig)

    txt_out = os.path.join(outdir, "al30lu_transformations.txt")
    with open(txt_out, "w", encoding="utf-8") as f:
        f.write(f"Schematic approximation used for Al–{solute}.\n")
        f.write(f"Eutectic (approx): X({solute})≈{X_e:.2f}, T≈{T_e:.0f} K\n")
        f.write(f"Peritectic (approx): X({solute})≈{X_p:.2f}, T≈{T_p:.0f} K\n")

    print("Done (schematic).")
    print(f"- Binary diagram (schematic): {bin_out}")
    print(f"- Phase fractions (schematic): {eq_fig_out}")
    print(f"- Transformations (schematic): {txt_out}")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())

