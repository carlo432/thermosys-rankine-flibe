#!/usr/bin/env python3
"""Moving-boundary boiler zones and FLiBe property uncertainty.

Run:  PYTHONPATH=. python3 runs/run_boiler_and_uncertainty.py
"""
import csv
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thermosys.analysis.uncertainty import DEFAULT_BANDS, corner_study, one_at_a_time
from thermosys.fluids import FlibeFluid, WaterFluid
from thermosys.transient import MovingBoundaryBoiler

warnings.filterwarnings("ignore")
OUT = pathlib.Path(__file__).parent
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
GRID = dict(color="#d8d8d4", linewidth=0.6)


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=11, loc="left", color="#0b0b0b")
    ax.set_xlabel(xlabel, fontsize=9, color="#52514e")
    ax.set_ylabel(ylabel, fontsize=9, color="#52514e")
    ax.grid(True, **GRID)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b8b8b4")
    ax.tick_params(colors="#52514e", labelsize=8)


def boiler_zone_map():
    f = WaterFluid()
    h_in = f.state_PQ(8e6, 0.0).h - 3.0e5
    qs = np.linspace(1.5e7, 2.5e8, 80)
    rows = []
    for q in qs:
        b = MovingBoundaryBoiler(fluid=f, P=8e6, mdot=50.0, h_in=h_in, q_total=float(q),
                                 length=20.0, M_total=500.0, n_cells=200)
        L1, L2, L3 = b.zone_lengths(b.steady_profile())
        rows.append((q, L1, L2, L3, b.outlet_quality(b.steady_profile())))

    with open(OUT / "boiler_zone_map.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["q_total_W", "L_subcooled_m", "L_twophase_m", "L_superheat_m", "x_out"])
        w.writerows(rows)

    q, L1, L2, L3, _ = map(np.array, zip(*rows))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.stackplot(q / 1e6, L1, L2, L3, colors=SERIES[:3], edgecolor="#fcfcfb", linewidth=0.8,
                 labels=["subcooled", "two-phase", "superheat"])
    style(ax, "Boiler zone lengths vs heat input (8 MPa, 50 kg/s)",
          "total heat input [MW]", "zone length [m]")
    ax.set_xlim(q[0] / 1e6, q[-1] / 1e6)
    ax.set_ylim(0, 20)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "boiler_zone_map.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    onset = q[np.argmax(L3 > 0)] / 1e6
    print(f"  boiler_zone_map: superheat onset at {onset:.1f} MW; "
          f"at {q[-1] / 1e6:.0f} MW zones are "
          f"{L1[-1]:.2f} / {L2[-1]:.2f} / {L3[-1]:.2f} m")


def boiler_transient():
    f = WaterFluid()
    h_in = f.state_PQ(8e6, 0.0).h - 3.0e5
    base = MovingBoundaryBoiler(fluid=f, P=8e6, mdot=50.0, h_in=h_in, q_total=1.4e8,
                                length=20.0, M_total=500.0, n_cells=200)
    steady = base.steady_profile()
    ramp = MovingBoundaryBoiler(fluid=f, P=8e6, mdot=50.0, h_in=h_in, q_total=1.4e8,
                                length=20.0, M_total=500.0, n_cells=200,
                                q_func=lambda t: 1.4e8 * (1.0 + 0.3 * (t > 5.0)))
    sol = ramp.solve(h0=steady, t_final=80.0)
    t = np.linspace(0, 80, 401)
    zones = np.array([ramp.zone_lengths(sol.sol(x)) for x in t])

    with open(OUT / "boiler_transient.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_s", "L_subcooled_m", "L_twophase_m", "L_superheat_m"])
        w.writerows(np.column_stack([t, zones]))

    fig, ax = plt.subplots(figsize=(7, 4))
    for i, (lab, c) in enumerate(zip(["subcooled", "two-phase", "superheat"], SERIES)):
        ax.plot(t, zones[:, i], color=c, lw=2, label=lab)
    ax.axvline(5.0, color="#52514e", lw=1, ls="--")
    ax.annotate("+30% heat input", xy=(5.0, 15), xytext=(8, 0),
                textcoords="offset points", fontsize=8, color="#52514e")
    style(ax, "Phase boundaries moving after a step in heat input", "time [s]",
          "zone length [m]")
    ax.legend(frameon=False, fontsize=9, loc="center right")
    fig.tight_layout()
    fig.savefig(OUT / "boiler_transient.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  boiler_transient: zones {zones[0].round(2)} -> {zones[-1].round(2)} m")


def uncertainty_figure():
    base = FlibeFluid()
    o = one_at_a_time(base)
    rows = [r for r in corner_study(base) if "T_return" in r]
    Ts = np.array([r["T_return"] for r in rows])

    with open(OUT / "flibe_uncertainty.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["property", "band", "T_return_low_K", "T_return_high_K", "range_K"])
        for n, v in o["effects"].items():
            w.writerow([n, DEFAULT_BANDS[n], v["low"], v["high"], v["range_K"]])

    names = sorted(o["effects"], key=lambda n: o["effects"][n]["range_K"])
    ranges = [o["effects"][n]["range_K"] for n in names]
    labels = [f"{n}  ±{DEFAULT_BANDS[n]:.0%}" for n in names]

    fig, ax = plt.subplots(figsize=(7, 3.6))
    bars = ax.barh(labels, ranges, color=SERIES[0], height=0.55)
    for b, r in zip(bars, ranges):
        ax.text(b.get_width() + 0.3, b.get_y() + b.get_height() / 2,
                f"{r:.2f} K", va="center", fontsize=8, color="#52514e")
    style(ax, "What moves the FLiBe return temperature",
          "swing in return temperature [K]", "")
    ax.text(0.0, 1.02, f"nominal {o['nominal']:.1f} K · all-corner spread "
                       f"{np.ptp(Ts):.1f} K · min freeze margin "
                       f"{min(r['freeze_margin'] for r in rows):.0f} K",
            transform=ax.transAxes, fontsize=8.5, color="#52514e")
    ax.set_xlim(0, max(ranges) * 1.3)
    fig.tight_layout()
    fig.savefig(OUT / "flibe_uncertainty.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  flibe_uncertainty: nominal {o['nominal']:.2f} K, corners "
          f"{Ts.min():.2f}..{Ts.max():.2f} K (spread {np.ptp(Ts):.2f} K)")
    for n in reversed(names):
        print(f"    {n:4} +/-{DEFAULT_BANDS[n]:.0%} -> {o['effects'][n]['range_K']:6.2f} K")


if __name__ == "__main__":
    print("boiler + uncertainty ->", OUT)
    boiler_zone_map()
    boiler_transient()
    uncertainty_figure()
