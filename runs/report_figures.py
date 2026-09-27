#!/usr/bin/env python3
"""Render the four report figures as compact SVG for embedding.

Run:  PYTHONPATH=. python3 runs/report_figures.py <outdir>
"""
import pathlib
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thermosys.analysis.uncertainty import DEFAULT_BANDS, corner_study, one_at_a_time
from thermosys.cycles import RankineCycle, RegenerativeCycle, ReheatCycle
from thermosys.fluids import FlibeFluid, WaterFluid
from thermosys.transient import FlibeLoopTransient, MovingBoundaryBoiler

warnings.filterwarnings("ignore")
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRIDC = "#14181c", "#5a636b", "#d9dcdf"
plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "font.family": "ui-monospace, Menlo, Consolas, monospace"})


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=10, loc="left", color=INK, pad=10)
    ax.set_xlabel(xlabel, fontsize=8.5, color=MUTED)
    ax.set_ylabel(ylabel, fontsize=8.5, color=MUTED)
    ax.grid(True, color=GRIDC, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b3b8bd")
    ax.tick_params(colors=MUTED, labelsize=8)


def save(fig, out, name):
    fig.tight_layout()
    fig.savefig(out / f"{name}.svg", format="svg", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}.svg  {(out / (name + '.svg')).stat().st_size // 1024} KB")


def fig_uncertainty(out):
    o = one_at_a_time(FlibeFluid())
    rows = [r for r in corner_study(FlibeFluid()) if "T_return" in r]
    Ts = np.array([r["T_return"] for r in rows])
    names = sorted(o["effects"], key=lambda n: o["effects"][n]["range_K"])
    ranges = [o["effects"][n]["range_K"] for n in names]
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    bars = ax.barh([f"{n}  ±{DEFAULT_BANDS[n]:.0%}" for n in names], ranges,
                   color=[SERIES[1] if r > 5 else "#9aa3aa" for r in ranges], height=0.5)
    for b, r in zip(bars, ranges):
        ax.text(b.get_width() + 0.35, b.get_y() + b.get_height() / 2, f"{r:.2f} K",
                va="center", fontsize=8, color=MUTED)
    style(ax, f"FLiBe return temperature · nominal {o['nominal']:.1f} K · "
              f"corner spread {np.ptp(Ts):.1f} K", "swing [K]", "")
    ax.set_xlim(0, max(ranges) * 1.32)
    save(fig, out, "fig_uncertainty")


def fig_zone_map(out):
    f = WaterFluid()
    h_in = f.state_PQ(8e6, 0.0).h - 3.0e5
    q = np.linspace(1.5e7, 2.5e8, 60)
    Z = []
    for qi in q:
        b = MovingBoundaryBoiler(fluid=f, P=8e6, mdot=50.0, h_in=h_in, q_total=float(qi),
                                 length=20.0, M_total=500.0, n_cells=200)
        Z.append(b.zone_lengths(b.steady_profile()))
    Z = np.array(Z)
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    ax.stackplot(q / 1e6, Z[:, 0], Z[:, 1], Z[:, 2], colors=SERIES[:3],
                 edgecolor="white", linewidth=0.7,
                 labels=["subcooled", "two-phase", "superheat"])
    style(ax, "Boiler phase boundaries vs heat input · 8 MPa, 50 kg/s",
          "heat input [MW]", "zone length [m]")
    ax.set_xlim(q[0] / 1e6, q[-1] / 1e6)
    ax.set_ylim(0, 20)
    ax.legend(frameon=False, fontsize=8, loc="upper right", labelcolor=MUTED)
    save(fig, out, "fig_zone_map")


def fig_reheat(out):
    Ps = np.geomspace(2e5, 6e6, 40)
    eta = np.array([ReheatCycle(P_reheat=float(p)).solve()["eta_th"] for p in Ps])
    plain = RankineCycle(P_boiler=8e6, P_cond=1e4, T_boiler=773.15).solve()["eta_th"]
    regen = RegenerativeCycle().solve()["eta_th"]
    i = int(eta.argmax())
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    ax.plot(Ps / 1e6, eta, color=SERIES[0], lw=2, label="reheat")
    ax.axhline(plain, color=SERIES[1], lw=1.6, ls="--", label="plain 8 MPa")
    ax.axhline(regen, color=SERIES[2], lw=1.6, ls=":", label="regenerative")
    ax.plot([Ps[i] / 1e6], [eta[i]], "o", ms=6, color=SERIES[0],
            markeredgecolor="white", markeredgewidth=1.8)
    ax.annotate(f"{Ps[i] / 1e6:.2f} MPa · {eta[i]:.4f}", xy=(Ps[i] / 1e6, eta[i]),
                xytext=(0, 10), textcoords="offset points", ha="center",
                fontsize=8, color=MUTED)
    ax.set_xscale("log")
    ax.set_xticks([0.2, 0.5, 1, 2, 4, 6])
    ax.set_xticklabels(["0.2", "0.5", "1", "2", "4", "6"])
    ax.minorticks_off()
    style(ax, "Reheat pressure optimum · 8 MPa / 500 °C / 10 kPa",
          "reheat pressure [MPa]", "thermal efficiency")
    ax.set_ylim(min(eta) - 0.002, max(regen, eta.max()) + 0.004)
    ax.legend(frameon=False, fontsize=8, loc="lower center", ncol=3, labelcolor=MUTED)
    save(fig, out, "fig_reheat")


def fig_loss_of_sink(out):
    net = FlibeLoopTransient(fluid=FlibeFluid(), mdot=3900.0, Q_core=500e6,
                             C_cold=13000000.0, T_sec_in=800.0, UA=4.5e6,
                             UA_func=lambda t: 0.0, deltaP=0.0)
    ref = FlibeLoopTransient(fluid=FlibeFluid(), mdot=3900.0, Q_core=500e6,
                             C_cold=13000000.0, T_sec_in=800.0, UA=4.5e6)
    sol = net.solve(ref.steady_guess(903.691595), 100.0)
    t = np.linspace(0, 100, 300)
    mean = np.array([(net.masses @ sol.sol(x)) / net.masses.sum() for x in t])
    slope = net.Q_core / (net.masses.sum() * net.cp)
    fig, ax = plt.subplots(figsize=(6.4, 3.1))
    ax.plot(t, mean, color=SERIES[0], lw=2.2, label="simulated mean")
    ax.plot(t, mean[0] + slope * t, color=SERIES[3], lw=1.6, ls="--",
            label=f"analytic Q/(M·cp) = {slope:.4f} K/s")
    style(ax, "Loss of heat sink · slope is an analytic consequence of the energy balance",
          "time [s]", "mass-weighted mean temperature [K]")
    ax.legend(frameon=False, fontsize=8, loc="upper left", labelcolor=MUTED)
    save(fig, out, "fig_loss_of_sink")


if __name__ == "__main__":
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out.mkdir(parents=True, exist_ok=True)
    print("report figures ->", out)
    fig_uncertainty(out)
    fig_zone_map(out)
    fig_reheat(out)
    fig_loss_of_sink(out)
