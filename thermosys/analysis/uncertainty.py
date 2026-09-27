"""Propagate FLiBe property uncertainty through the loop.

Verification says the equations are solved right. It says nothing about whether the
property correlations feeding them are right. For FLiBe that gap is real: the viscosity
fit comes from a 873-1073 K dataset and is used from 732 to 1703 K, and the conductivity
line is fitted over 800-1200 K. This module answers the only question that matters for
using the results -- how far does the answer move when the properties move?

Bands used (literature scatter for FLiBe, deliberately generous rather than optimistic):
    mu  +/-20%      the widest-scattering property between sources
    k   +/-15%
    rho +/-2%       well constrained
    cp  +/-3%
"""
import itertools

from thermosys.components.htc import loop_UA
from thermosys.cycles import FlibeLoop
from thermosys.fluids.perturbed import PerturbedFluid

DEFAULT_BANDS = dict(mu=0.20, k=0.15, rho=0.02, cp=0.03)


def solve_with_computed_UA(fluid, mdot=3900.0, Q_core=500e6, secondary=None,
                           mdot_sec=2500.0, T_sec_in=800.0, T_ref=900.0,
                           D=0.02, area=900.0, n_tubes=2000, deltaP=1.0e5):
    """Solve the loop with UA computed from the coolant's own properties.

    UA is evaluated once at T_ref rather than at the converged return temperature: the
    loop would otherwise need an inner property iteration, and the resulting UA shift is
    small next to the property uncertainty this module is quantifying.
    """
    from thermosys.fluids import ConstantCpFluid
    secondary = secondary or ConstantCpFluid(cp=5200.0)
    UA = loop_UA(fluid, T_ref, 1.0e5, mdot, D, area, n_tubes=n_tubes)
    loop = FlibeLoop(fluid, mdot, Q_core, secondary, mdot_sec, T_sec_in, UA,
                     deltaP=deltaP)
    result = loop.solve()
    result["UA"] = UA
    return result


def corner_study(base_fluid, bands=None, **kw):
    """Evaluate every combination of property extremes (2^4 = 16 corners) plus nominal.

    Corners rather than Monte Carlo: with four inputs and a monotone response the extreme
    combinations bound the range exactly, and the result is reproducible without a seed.
    """
    bands = bands or DEFAULT_BANDS
    names = sorted(bands)
    rows = []
    for signs in itertools.product((-1, 0, 1), repeat=len(names)):
        if 0 in signs and any(s != 0 for s in signs):
            continue                          # nominal only as the all-zero case
        factors = {n: 1.0 + s * bands[n] for n, s in zip(names, signs)}
        fluid = PerturbedFluid(base_fluid, **factors)
        try:
            r = solve_with_computed_UA(fluid, **kw)
        except Exception as exc:              # a corner may leave the liquid range
            rows.append(dict(factors=factors, error=str(exc)[:80]))
            continue
        rows.append(dict(factors=factors, T_return=r["T_return"],
                         T_core_out=r["T_core_out"], UA=r["UA"], eps=r["eps"],
                         freeze_margin=r["freeze_margin"]))
    return rows


def one_at_a_time(base_fluid, bands=None, **kw):
    """Sensitivity of T_return to each property alone: the ranking that tells you which
    correlation is worth improving."""
    bands = bands or DEFAULT_BANDS
    nominal = solve_with_computed_UA(PerturbedFluid(base_fluid), **kw)
    out = {"nominal": nominal["T_return"], "nominal_UA": nominal["UA"], "effects": {}}
    for name, band in sorted(bands.items()):
        swing = {}
        for sign, label in ((-1, "low"), (1, "high")):
            fluid = PerturbedFluid(base_fluid, **{name: 1.0 + sign * band})
            swing[label] = solve_with_computed_UA(fluid, **kw)["T_return"]
        swing["range_K"] = abs(swing["high"] - swing["low"])
        out["effects"][name] = swing
    return out
