"""Lumped control-volume transients.

Every right-hand side here is a conservation equation whose terms depend on the state.
That is the whole point: a model that relaxes onto a precomputed steady answer has no
test of them could fail.  Deleting any term below breaks a named test -- see
FALSIFIABILITY.md.

Scope is deliberately bounded: lumped CVs, quasi-static turbomachinery, no
moving-boundary or drum-level two-phase boiler.
"""
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from thermosys.components.hx import hx_duty


def integrate(rhs, t_span, y0, **kwargs):
    """Stiff integration. BDF is not a preference: these systems span time constants from
    a per-CV residence time to a full loop transit."""
    opts = dict(method="BDF", rtol=1e-10, atol=1e-8, dense_output=True)
    opts.update(kwargs)
    return solve_ivp(rhs, t_span, np.asarray(y0, dtype=float), **opts)


# ---------------------------------------------------------------------------
# 6.1 single control volume, analytic anchor
# ---------------------------------------------------------------------------
@dataclass
class SingleCV:
    """M*cp*dT/dt = mdot*cp*(T_in - T) + Q."""
    M: float
    cp: float
    mdot: float
    T_in: float
    Q: float

    @property
    def tau(self):
        return self.M / self.mdot

    @property
    def T_ss(self):
        return self.T_in + self.Q / (self.mdot * self.cp)

    def rhs(self, t, y):
        return [(self.mdot * self.cp * (self.T_in - y[0]) + self.Q) / (self.M * self.cp)]

    def analytic(self, t, T0):
        return self.T_ss + (T0 - self.T_ss) * math.exp(-t / self.tau)

    def solve(self, T0, t_final, **kwargs):
        return integrate(self.rhs, (0.0, t_final), [T0], **kwargs)


# ---------------------------------------------------------------------------
# 6.2 transport delay: N CVs in series
# ---------------------------------------------------------------------------
def erlang_response(t, N, tau_cv):
    """Outlet response of N well-mixed CVs in series to a unit inlet step."""
    z = t / tau_cv
    return 1.0 - math.exp(-z) * sum(z ** k / math.factorial(k) for k in range(N))


@dataclass
class CVChain:
    """N equal CVs in series; inlet steps from 0 to 1 at t = 0. dT_i/dt = (T_{i-1}-T_i)/tau."""
    N: int
    tau_cv: float

    def rhs(self, t, y):
        upstream = np.concatenate(([1.0], y[:-1]))
        return (upstream - y) / self.tau_cv

    def solve(self, t_final, **kwargs):
        return integrate(self.rhs, (0.0, t_final), np.zeros(self.N), **kwargs)


# ---------------------------------------------------------------------------
# 6.3 FLiBe loop network
# ---------------------------------------------------------------------------
@dataclass
class FlibeLoopTransient:
    """Node chain: core -> hot leg (n_leg CVs) -> HX -> cold leg (n_leg CVs) -> core.

    Energy balance per node, constant cp:
        M_i*cp*dT_i/dt = mdot*cp*(T_upstream - T_i) + S_i
    with S_core = Q_core + W_pump and S_hx = -Q_hx.  Q_hx uses the exchanger INLET
    temperature, which is what makes this network's steady state identical to the
    algebraic FlibeLoop solution -- the steady-state test depends on that.
    """
    fluid: object
    mdot: float
    Q_core: float
    C_cold: float
    T_sec_in: float
    UA: float
    M_core: float = 20000.0
    M_leg: float = 10000.0
    M_hx: float = 20000.0
    n_leg: int = 5
    deltaP: float = 1.0e5
    P: float = 1.0e5
    mdot_func: object = None
    UA_func: object = None
    Q_func: object = None

    masses: np.ndarray = field(init=False)
    i_core: int = field(init=False)
    i_hx: int = field(init=False)

    def __post_init__(self):
        self.i_core = 0
        self.i_hx = 1 + self.n_leg
        self.masses = np.array(
            [self.M_core] + [self.M_leg] * self.n_leg + [self.M_hx]
            + [self.M_leg] * self.n_leg, dtype=float)

    @property
    def n_nodes(self):
        return len(self.masses)

    @property
    def cp(self):
        # Evaluated at the coolant's own reference temperature, not the secondary inlet:
        # the secondary can legitimately sit outside the FLiBe liquid range.
        return self.fluid.cp(self.fluid.T_min, self.P)

    def _rho(self, T):
        """Properties are evaluated on a clamped temperature so the integrator may probe
        outside the liquid range without crashing; excursions are caught on the solution."""
        return self.fluid.rho(min(max(T, self.fluid.T_min), self.fluid.T_max), self.P)

    def mdot_at(self, t):
        return self.mdot if self.mdot_func is None else self.mdot_func(t)

    def UA_at(self, t):
        return self.UA if self.UA_func is None else self.UA_func(t)

    def Q_at(self, t):
        return self.Q_core if self.Q_func is None else self.Q_func(t)

    def pump_work(self, T_return, mdot):
        return mdot * self.deltaP / self._rho(T_return)

    def rhs(self, t, y):
        cp = self.cp
        mdot = self.mdot_at(t)
        C_h = mdot * cp
        upstream = np.roll(y, 1)          # closed loop: node 0's upstream is the last node
        dTdt = mdot * cp * (upstream - y)

        dTdt[self.i_core] += self.Q_at(t) + self.pump_work(y[-1], mdot)

        UA = self.UA_at(t)
        if UA > 0.0 and C_h > 0.0:
            q = hx_duty(y[self.i_hx - 1], self.T_sec_in, C_h, self.C_cold, UA)
            dTdt[self.i_hx] -= q
        return dTdt / (self.masses * cp)

    def steady_guess(self, T_return):
        cp = self.cp
        rise = (self.Q_core + self.pump_work(T_return, self.mdot)) / (self.mdot * cp)
        return np.array([T_return + rise] * (1 + self.n_leg) + [T_return] * (1 + self.n_leg))

    def solve(self, y0, t_final, check_bounds=True, **kwargs):
        sol = integrate(self.rhs, (0.0, t_final), y0, **kwargs)
        if check_bounds:
            lo, hi = float(np.min(sol.y)), float(np.max(sol.y))
            if lo < self.fluid.T_min:
                raise ValueError(
                    f"FLiBe froze during the transient: minimum node temperature {lo:.2f} K "
                    f"is below the {self.fluid.T_min:.1f} K freezing point")
            if hi > self.fluid.T_max:
                raise ValueError(
                    f"FLiBe boiled during the transient: maximum node temperature {hi:.2f} K "
                    f"is above the {self.fluid.T_max:.1f} K boiling point")
        sol.freeze_margin = float(np.min(sol.y) - self.fluid.T_min)
        sol.T_core_out = sol.y[self.i_core]
        sol.T_return = sol.y[-1]
        sol.mean_T = (self.masses @ sol.y) / self.masses.sum()
        return sol
