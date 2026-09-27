#!/usr/bin/env python3
"""Transient scenarios: CSV + figure per scenario, written next to this script.

Run:  python3 runs/run_transient_scenarios.py
"""
import csv
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thermosys.cycles import FlibeLoop
from thermosys.fluids import ConstantCpFluid, FlibeFluid
from thermosys.transient import FlibeLoopTransient

warnings.filterwarnings("ignore")
OUT = pathlib.Path(__file__).parent

# categorical slots 1-4, validated for CVD separation on a light surface
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
GRID = dict(color="#d8d8d4", linewidth=0.6)


def network(**kw):
    base = dict(fluid=FlibeFluid(), mdot=3900.0, Q_core=500e6,
                C_cold=13000000.0, T_sec_in=800.0, UA=4.5e6)
    return FlibeLoopTransient(**{**base, **kw})


def algebraic(Q=500e6, mdot=3900.0, UA=4.5e6, deltaP=1.0e5):
    return FlibeLoop(FlibeFluid(), mdot, Q, ConstantCpFluid(cp=5200.0), 2500.0, 800.0,
                     UA, deltaP=deltaP).solve()


def style(ax, title, xlabel="time [s]", ylabel="temperature [K]"):
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


def write(name, header, rows):
    with open(OUT / f"{name}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"  wrote {name}.csv / {name}.png")


def scenario_power_step():
    name = "scenario_power_step"
    net = network(Q_core=550e6)
    sol = net.solve(network().steady_guess(903.691595), 4000.0)
    ref = algebraic(Q=550e6)
    t = np.linspace(0, 4000, 601)
    core = np.array([sol.sol(x)[net.i_core] for x in t])
    ret = np.array([sol.sol(x)[-1] for x in t])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, core, color=SERIES[0], lw=2, label="core outlet")
    ax.plot(t, ret, color=SERIES[1], lw=2, label="loop return")
    ax.axhline(ref["T_return"], color="#52514e", lw=1, ls="--")
    ax.annotate(f"algebraic solve  {ref['T_return']:.2f} K",
                xy=(t[-1], ref["T_return"]), xytext=(-8, 8),
                textcoords="offset points", ha="right", fontsize=8, color="#52514e")
    style(ax, "Step to 550 MW core power — settles on the steady-state solve")
    ax.legend(frameon=False, fontsize=9, loc="center right")
    save(fig, name)
    write(name, ["t_s", "T_core_out_K", "T_return_K"], zip(t, core, ret))
    return f"final return {ret[-1]:.4f} K vs algebraic {ref['T_return']:.4f} K"


def scenario_flow_step():
    name = "scenario_flow_step"
    hi, lo = network(deltaP=0.0), network(mdot=3120.0, deltaP=0.0)
    r_hi = hi.solve(hi.steady_guess(900.0), 6000.0)
    r_lo = lo.solve(lo.steady_guess(900.0), 6000.0)
    t = np.linspace(0, 6000, 601)
    rise_hi = np.array([r_hi.sol(x)[hi.i_core] - r_hi.sol(x)[-1] for x in t])
    rise_lo = np.array([r_lo.sol(x)[lo.i_core] - r_lo.sol(x)[-1] for x in t])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, rise_hi, color=SERIES[0], lw=2, label="mdot = 3900 kg/s")
    ax.plot(t, rise_lo, color=SERIES[1], lw=2, label="mdot = 3120 kg/s (-20%)")
    style(ax, "Core temperature rise scales as 1/mdot", ylabel="core rise [K]")
    ax.legend(frameon=False, fontsize=9)
    save(fig, name)
    write(name, ["t_s", "rise_3900_K", "rise_3120_K"], zip(t, rise_hi, rise_lo))
    return f"rise ratio {rise_lo[-1] / rise_hi[-1]:.6f} vs exact {3900 / 3120:.6f}"


def scenario_coastdown():
    name = "scenario_coastdown"
    net = network(mdot_func=lambda x: 3900.0 * max(np.exp(-x / 5.0), 0.05))
    sol = net.solve(network().steady_guess(903.691595), 20.0)
    t = np.linspace(0, 20, 401)
    core = np.array([sol.sol(x)[net.i_core] for x in t])
    mdot = np.array([net.mdot_at(x) for x in t])
    peak, t_peak = float(core.max()), float(t[core.argmax()])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    ax1.plot(t, core, color=SERIES[0], lw=2)
    ax1.plot([t_peak], [peak], "o", ms=7, color=SERIES[0],
             markeredgecolor="#fcfcfb", markeredgewidth=2)
    # still climbing at t_final: this is the maximum in the window, not a turning point
    ax1.annotate(f"max in window {peak:.1f} K", xy=(t_peak, peak),
                 xytext=(-12, -16), textcoords="offset points", ha="right",
                 fontsize=8, color="#52514e")
    style(ax1, "Pump coastdown: core outlet rises as flow decays", xlabel="")
    ax2.plot(t, mdot, color=SERIES[1], lw=2)
    style(ax2, "", ylabel="mdot [kg/s]")
    save(fig, name)
    write(name, ["t_s", "T_core_out_K", "mdot_kg_s"], zip(t, core, mdot))
    return (f"max core outlet {peak:.4f} K at t = {t_peak:.3f} s "
            f"(still rising at the end of the window)")


def scenario_loss_of_sink():
    name = "scenario_loss_of_sink"
    net = network(UA_func=lambda x: 0.0, deltaP=0.0)
    sol = net.solve(network().steady_guess(903.691595), 100.0)
    t = np.linspace(0, 100, 401)
    mean = np.array([(net.masses @ sol.sol(x)) / net.masses.sum() for x in t])
    expected = net.Q_core / (net.masses.sum() * net.cp)
    slope = (mean[-1] - mean[0]) / (t[-1] - t[0])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, mean, color=SERIES[0], lw=2, label="mass-weighted mean")
    ax.plot(t, mean[0] + expected * t, color=SERIES[3], lw=2, ls="--",
            label=f"analytic Q/(M cp) = {expected:.4f} K/s")
    style(ax, "Loss of heat sink — adiabatic slope is an analytic consequence")
    ax.legend(frameon=False, fontsize=9)
    save(fig, name)
    write(name, ["t_s", "mean_T_K"], zip(t, mean))
    return f"slope {slope:.6f} K/s vs analytic {expected:.6f} K/s"


if __name__ == "__main__":
    print("transient scenarios ->", OUT)
    for fn in (scenario_power_step, scenario_flow_step, scenario_coastdown,
               scenario_loss_of_sink):
        print(f"  {fn.__name__}: {fn()}")
