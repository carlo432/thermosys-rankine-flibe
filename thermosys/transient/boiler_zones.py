"""Once-through boiler with moving phase boundaries.

WHAT THIS IS, precisely, because the name is easy to over-claim:

  A distributed enthalpy model of the tube -- N finite volumes, each with a real energy
  balance -- on top of which the subcooled / two-phase / superheat boundaries are
  *diagnosed* by locating where the solved enthalpy field crosses h_f and h_g, with
  linear interpolation inside the crossing cell.  The boundaries move because the
  enthalpy field moves.

WHAT THIS IS NOT:

  The classical three-state moving-boundary formulation, which carries zone lengths as
  ODE states.  That form goes singular when a zone vanishes (its length ODE divides by
  the zone's own length), which is the standard reason those models stall.  Diagnosing
  boundaries from a distributed field costs more states but degrades gracefully: a zone
  that vanishes simply reports zero length, and the test suite exercises exactly that.

BOUNDING SIMPLIFICATION, stated rather than buried: cell inventory is uniform and
constant (M_i = M_total/N).  Real vapour cells hold far less mass than liquid cells, and
capturing that requires a coupled mass balance that redistributes flow along the tube --
genuinely out of scope here.  The consequence is that transport delay is uniform along
the tube rather than accelerating through the two-phase region.  Steady-state zone
lengths are unaffected, which is why they can be checked against a closed form.
"""
from dataclasses import dataclass, field

import numpy as np

from thermosys.transient.lumped import integrate


@dataclass
class MovingBoundaryBoiler:
    """Once-through tube at fixed pressure with a specified heat input per unit length."""
    fluid: object
    P: float
    mdot: float
    h_in: float
    q_total: float          # W over the whole tube
    length: float = 20.0
    M_total: float = 500.0
    n_cells: int = 60
    q_func: object = None    # optional q_total(t)
    mdot_func: object = None

    h_f: float = field(init=False)
    h_g: float = field(init=False)

    def __post_init__(self):
        self.h_f = self.fluid.state_PQ(self.P, 0.0).h
        self.h_g = self.fluid.state_PQ(self.P, 1.0).h

    # -- inputs ---------------------------------------------------------------
    def q_at(self, t):
        return self.q_total if self.q_func is None else self.q_func(t)

    def mdot_at(self, t):
        return self.mdot if self.mdot_func is None else self.mdot_func(t)

    @property
    def cell_mass(self):
        return self.M_total / self.n_cells

    @property
    def dz(self):
        return self.length / self.n_cells

    # -- dynamics -------------------------------------------------------------
    def rhs(self, t, h):
        """M_i dh_i/dt = mdot (h_{i-1} - h_i) + q_cell, with h_0 the tube inlet."""
        mdot = self.mdot_at(t)
        q_cell = self.q_at(t) / self.n_cells
        upstream = np.concatenate(([self.h_in], h[:-1]))
        return (mdot * (upstream - h) + q_cell) / self.cell_mass

    # -- steady state, in closed form -----------------------------------------
    def steady_profile(self):
        """Steady cell values: mdot(h_{i-1} - h_i) + q_cell = 0 marched from the inlet."""
        q_cell = self.q_total / self.n_cells
        return self.h_in + q_cell * (np.arange(self.n_cells) + 1.0) / self.mdot

    def analytic_zone_lengths(self):
        """L1 = mdot(h_f - h_in)/q', L2 = mdot(h_g - h_f)/q', clipped to the tube."""
        q_prime = self.q_total / self.length
        L1 = self.mdot * (self.h_f - self.h_in) / q_prime
        L2 = self.mdot * (self.h_g - self.h_f) / q_prime
        L1 = min(max(L1, 0.0), self.length)
        L2 = min(max(L2, 0.0), self.length - L1)
        return L1, L2, self.length - L1 - L2

    # -- boundary diagnosis ---------------------------------------------------
    def _crossing(self, h, target):
        """Axial position where the enthalpy profile crosses `target`.

        In upwind finite volumes h_i is both the cell average and the value leaving the
        cell, so it belongs at the DOWNSTREAM face z = (i+1)*dz; the inlet face carries
        h_in.  Interpolating on that grid makes the diagnosis exact for the linear steady
        profile, independent of cell count.

        Returns 0.0 if the tube is already past `target` at inlet, and `length` if it
        never gets there -- so a vanished zone reports zero length rather than dividing
        by it, which is where the classical three-state formulation goes singular.
        """
        z = np.arange(self.n_cells + 1) * self.dz          # inlet face .. outlet face
        hf = np.concatenate(([self.h_in], np.asarray(h)))
        if hf[0] >= target:
            return 0.0
        if hf[-1] < target:
            return self.length
        i = int(np.argmax(hf >= target))
        h_lo, h_hi = hf[i - 1], hf[i]
        if h_hi == h_lo:
            return float(z[i])
        return float(z[i - 1] + (target - h_lo) / (h_hi - h_lo) * self.dz)

    def zone_lengths(self, h):
        """(subcooled, two-phase, superheat) lengths diagnosed from an enthalpy field."""
        z_f = self._crossing(np.asarray(h), self.h_f)
        z_g = self._crossing(np.asarray(h), self.h_g)
        z_g = max(z_g, z_f)
        return z_f, z_g - z_f, self.length - z_g

    def outlet_quality(self, h):
        h_out = float(np.asarray(h)[-1])
        if h_out <= self.h_f:
            return 0.0
        if h_out >= self.h_g:
            return 1.0
        return (h_out - self.h_f) / (self.h_g - self.h_f)

    # -- solve ----------------------------------------------------------------
    def solve(self, h0=None, t_final=400.0, **kwargs):
        if h0 is None:
            h0 = np.full(self.n_cells, self.h_in)
        sol = integrate(self.rhs, (0.0, t_final), h0, **kwargs)
        sol.zone_lengths = np.array([self.zone_lengths(sol.y[:, k])
                                     for k in range(sol.y.shape[1])])
        sol.h_out = sol.y[-1]
        sol.x_out = np.array([self.outlet_quality(sol.y[:, k])
                              for k in range(sol.y.shape[1])])
        return sol
