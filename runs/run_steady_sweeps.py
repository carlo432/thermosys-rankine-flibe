#!/usr/bin/env python3
"""Steady-state sweeps: FLiBe loop UA sensitivity and reheat-pressure optimum.

Run:  python3 runs/run_steady_sweeps.py
"""
import csv
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thermosys.cycles import FlibeLoop, RankineCycle, ReheatCycle, RegenerativeCycle
from thermosys.fluids import ConstantCpFluid, FlibeFluid

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
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b8b8b4")
    ax.tick_params(colors="#52514e", labelsize=8)


def loop_ua_sweep():
    UAs = np.linspace(2.0e6, 1.2e7, 60)
    rows = []
    for UA in UAs:
        r = FlibeLoop(FlibeFluid(), 3900.0, 500e6, ConstantCpFluid(cp=5200.0),
                      2500.0, 800.0, float(UA)).solve()
        rows.append((UA, r["T_return"], r["T_core_out"], r["eps"]))
    with open(OUT / "sweep_loop_UA.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["UA_W_K", "T_return_K", "T_core_out_K", "effectiveness"])
        w.writerows(rows)

    UA, Tret, Thot, eps = map(np.array, zip(*rows))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(UA / 1e6, Tret, color=SERIES[0], lw=2, label="loop return")
    ax.plot(UA / 1e6, Thot, color=SERIES[1], lw=2, label="core outlet")
    ax.axhline(732.0, color="#e34948", lw=1.5, ls="--")
    ax.annotate("FLiBe freezing 732 K", xy=(UA[-1] / 1e6, 732.0), xytext=(-6, 6),
                textcoords="offset points", ha="right", fontsize=8, color="#e34948")
    style(ax, "FLiBe loop temperatures vs exchanger size",
          "UA [MW/K]", "temperature [K]")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "sweep_loop_UA.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  sweep_loop_UA: T_return {Tret[0]:.2f} -> {Tret[-1]:.2f} K over "
          f"UA {UA[0]:.2g} -> {UA[-1]:.2g} W/K")


def reheat_pressure_sweep():
    Ps = np.geomspace(2.0e5, 6.0e6, 50)
    etas, x6s = [], []
    for P in Ps:
        r = ReheatCycle(P_reheat=float(P)).solve()
        etas.append(r["eta_th"])
        x6s.append(r["x6"])
    etas, x6s = np.array(etas), np.array(x6s)
    plain = RankineCycle(P_boiler=8e6, P_cond=1e4, T_boiler=773.15).solve()["eta_th"]
    regen = RegenerativeCycle().solve()["eta_th"]
    best = int(etas.argmax())

    with open(OUT / "sweep_reheat_pressure.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["P_reheat_Pa", "eta_th", "x6"])
        w.writerows(zip(Ps, etas, x6s))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(Ps / 1e6, etas, color=SERIES[0], lw=2, label="reheat cycle")
    ax.axhline(plain, color=SERIES[1], lw=2, ls="--", label="plain 8 MPa cycle")
    ax.axhline(regen, color=SERIES[2], lw=2, ls=":", label="regenerative (1 open FWH)")
    ax.plot([Ps[best] / 1e6], [etas[best]], "o", ms=7, color=SERIES[0],
            markeredgecolor="#fcfcfb", markeredgewidth=2)
    ax.annotate(f"optimum {Ps[best] / 1e6:.2f} MPa, {etas[best]:.4f}",
                xy=(Ps[best] / 1e6, etas[best]), xytext=(0, 12),
                textcoords="offset points", ha="center", fontsize=8, color="#52514e")
    ax.set_xscale("log")
    ticks = [0.2, 0.5, 1.0, 2.0, 4.0, 6.0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(top=max(etas.max(), regen) + 0.004)
    style(ax, "Reheat pressure optimum, 8 MPa / 500 C / 10 kPa",
          "reheat pressure [MPa]", "thermal efficiency")
    ax.legend(frameon=False, fontsize=9, loc="lower center")
    fig.tight_layout()
    fig.savefig(OUT / "sweep_reheat_pressure.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  sweep_reheat_pressure: optimum {Ps[best]:.4g} Pa at eta {etas[best]:.6f} "
          f"(plain {plain:.6f}, regenerative {regen:.6f})")


if __name__ == "__main__":
    print("steady sweeps ->", OUT)
    loop_ua_sweep()
    reheat_pressure_sweep()
