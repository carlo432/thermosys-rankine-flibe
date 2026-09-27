"""Counterflow heat exchanger, effectiveness-NTU."""
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class HXResult:
    hot_out: object
    cold_out: object
    q: float
    eps: float
    NTU: float
    Cr: float
    C_hot: float
    C_cold: float

    def __iter__(self):
        """Backwards-compatible (hot_out, cold_out, q) unpacking."""
        return iter((self.hot_out, self.cold_out, self.q))


def effectiveness_counterflow(NTU, Cr):
    """Counterflow effectiveness. Cr -> 1 is a removable singularity, handled explicitly."""
    if NTU < 0.0:
        raise ValueError(f"NTU must be non-negative, got {NTU}")
    if not 0.0 <= Cr <= 1.0:
        raise ValueError(f"Cr must lie in [0, 1], got {Cr}")
    if abs(Cr - 1.0) < 1e-12:
        return NTU / (1.0 + NTU)
    e = math.exp(-NTU * (1.0 - Cr))
    return (1.0 - e) / (1.0 - Cr * e)


def counterflow_hx(hot_in, cold_in, mdot_h, mdot_c, UA):
    """Solve a counterflow exchanger. Returns HXResult (unpacks as hot_out, cold_out, q)."""
    C_h = mdot_h * hot_in.fluid.cp(hot_in.T, hot_in.P)
    C_c = mdot_c * cold_in.fluid.cp(cold_in.T, cold_in.P)
    C_min, C_max = min(C_h, C_c), max(C_h, C_c)
    Cr = C_min / C_max
    NTU = UA / C_min
    eps = effectiveness_counterflow(NTU, Cr)
    q = eps * C_min * (hot_in.T - cold_in.T)
    hot_out = hot_in.fluid.state_TP(hot_in.T - q / C_h, hot_in.P)
    cold_out = cold_in.fluid.state_TP(cold_in.T + q / C_c, cold_in.P)
    return HXResult(hot_out, cold_out, q, eps, NTU, Cr, C_h, C_c)


def hx_duty(T_hot_in, T_cold_in, C_h, C_c, UA):
    """Duty only, from temperatures and capacity rates. Used inside ODE right-hand sides
    where allocating State objects every step is wasteful."""
    C_min, C_max = min(C_h, C_c), max(C_h, C_c)
    eps = effectiveness_counterflow(UA / C_min, C_min / C_max)
    return eps * C_min * (T_hot_in - T_cold_in)
