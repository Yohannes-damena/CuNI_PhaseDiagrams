import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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


def stable_phase_set_at(dbf: Database, comps: list[str], phases: list[str], x_lu: float, temp_k: float,
                        threshold: float = 1e-6) -> set[str]:
    eqt = equilibrium(dbf, comps, phases, {v.X("LU"): x_lu, v.T: temp_k, v.P: 101325})
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


def classify_invariant(dbf: Database, comps: list[str], phases: list[str], x_lu: float, temp_k: float,
                       delta_k: float = 2.0) -> str:
    above = stable_phase_set_at(dbf, comps, phases, x_lu, temp_k + delta_k)
    below = stable_phase_set_at(dbf, comps, phases, x_lu, temp_k - delta_k)
    if ("LIQUID" in above and len(above) == 2) and ("LIQUID" not in below and len(below) >= 2):
        return "eutectic"
    if ("LIQUID" in above and len(above) == 2) and ("LIQUID" not in below and len(below) == 1):
        return "peritectic"
    return "invariant"


def plot_binary_with_annotations(dbf: Database, components: list[str], phases: list[str],
                                 tmin: float, tmax: float, outpath: str) -> None:
    binplot(
        dbf,
        components,
        phases,
        conditions={
            v.X("LU"): np.linspace(0.0, 1.0, 141),
            v.T: np.linspace(tmin, tmax, 141),
            v.P: 101325,
        },
    )
    ax = plt.gca()
    ax.set_xlim(0, 1)
    ax.set_xlabel("X(LU)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Al–Lu Binary Phase Diagram")

    xs = np.linspace(0, 1, 81)
    ts = np.linspace(tmin, tmax, 161)
    eq = equilibrium(dbf, components, phases, {v.X("LU"): xs, v.T: ts, v.P: 101325})
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
            itype = classify_invariant(dbf, components, phases, Xc, Tc)
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
                          x_target: float, tmin: float, tmax: float,
                          fig_outpath: str, txt_outpath: str) -> None:
    T_range = (tmin, tmax, 161)
    EQ = equilibrium(dbf, components, phases, {v.X("LU"): x_target, v.T: T_range, v.P: 101325})
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
    ax.set_title(f"Al-{int(round(100*x_target))}Lu Phase Fractions vs Temperature")
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
        itype = classify_invariant(dbf, components, phases, x_target, Tk)
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
            f.write(f"Detected invariant transformations for X(LU)={x_target:.2f}:\n")
            for itype, Tk in annotated:
                f.write(f"  - {itype.capitalize()} near T ≈ {Tk:.1f} K ({Tk - 273.15:.1f} °C)\n")
        else:
            f.write("No clear invariant transformation detected along this exact composition path.\n")


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    db_path = args.db
    if not os.path.exists(db_path):
        print(f"Error: database not found at '{db_path}'", file=sys.stderr)
        return 1

    dbf = Database(db_path)
    components = ["AL", "LU", "VA"]
    phases = get_phase_list(dbf, args.phases)

    # Informational note for SGTE unary databases
    print("Note: If your TDB is a unary (pure-elements) database, the diagram may lack Al–Lu intermetallics and "
          "true eutectic/peritectic features. Use an assessed Al–Lu database for realistic results.\n")

    # Binary diagram with annotations
    binary_out = os.path.join(args.outdir, "al_lu_binary.png")
    print(f"Generating binary phase diagram → {binary_out}")
    plot_binary_with_annotations(dbf, components, phases, args.binplot_tmin, args.binplot_tmax, binary_out)

    # Equilibrium path at fixed composition
    eq_fig_out = os.path.join(args.outdir, "al30lu_phasefractions.png")
    eq_txt_out = os.path.join(args.outdir, "al30lu_transformations.txt")
    print(f"Computing Al-{int(round(100*args.xlu))}Lu equilibrium and phase fractions → {eq_fig_out}")
    plot_equilibrium_path(dbf, components, phases, args.xlu, args.tmin, args.tmax, eq_fig_out, eq_txt_out)

    print("Done.")
    print(f"- Binary diagram: {binary_out}")
    print(f"- Phase fractions: {eq_fig_out}")
    print(f"- Transformations summary: {eq_txt_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

