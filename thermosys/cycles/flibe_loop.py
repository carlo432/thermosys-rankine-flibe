"""Closed single-phase FLiBe loop: pump -> core -> counterflow HX -> return.

The loop closes on itself, so the return temperature is the root of a
residual, found with brentq over a scanned bracket.  Q_hx comes from the heat exchanger
model and is NOT defined as Q_core - W_pump (rule R2) -- the energy balance is therefore
a real check on the solution rather than a restatement of it.
"""
from scipy.optimize import brentq

from thermosys.components import counterflow_hx, heat_source
from thermosys.components.hx import hx_duty


class LoopDomainError(ValueError):
    """No loop solution exists inside the coolant's liquid range."""


class FlibeLoop:
    def __init__(self, fluid, mdot, Q_core, secondary, mdot_sec, T_sec_in, UA,
                 deltaP=1.0e5, P=1.0e5):
        self.fluid = fluid
        self.mdot = mdot
        self.Q_core = Q_core
        self.secondary = secondary
        self.mdot_sec = mdot_sec
        self.T_sec_in = T_sec_in
        self.UA = UA
        self.deltaP = deltaP
        self.P = P

    # -- physics ---------------------------------------------------------------
    def pump_work(self, T_return):
        """Incompressible pump work [W]. FLiBe has no EOS, so this is dP/rho, not dh."""
        return self.mdot * self.deltaP / self.fluid.rho(T_return, self.P)

    def _core_outlet(self, T_return):
        cold = self.fluid.state_TP(T_return, self.P)
        return heat_source(cold, self.mdot, self.Q_core + self.pump_work(T_return))

    def _exchanger(self, T_return):
        hot_in = self._core_outlet(T_return)
        cold_in = self.secondary.state_TP(self.T_sec_in, self.P)
        return counterflow_hx(hot_in, cold_in, self.mdot, self.mdot_sec, self.UA)

    def residual(self, T_return):
        """Mismatch between the HX hot-side outlet and the assumed return temperature.

        Evaluated on clamped properties and bare temperatures rather than State objects:
        while bracketing, trial temperatures legitimately fall outside the liquid range and
        must return a finite residual instead of raising, or the bracket scan cannot tell
        which side of the range the root lies on.  The duty still comes from the exchanger
        model (rule R2) -- `hx_duty` is the same effectiveness relation `counterflow_hx`
        uses, and `counterflow_hx` itself produces the final states in `solve`.
        """
        cp = self.fluid.cp(self._clamp(T_return), self.P)
        C_h = self.mdot * cp
        C_c = self.mdot_sec * self.secondary.cp(self.T_sec_in, self.P)
        W_pump = self.mdot * self.deltaP / self.fluid.rho(self._clamp(T_return), self.P)
        T_hot = T_return + (self.Q_core + W_pump) / C_h
        q = hx_duty(T_hot, self.T_sec_in, C_h, C_c, self.UA)
        return T_hot - q / C_h - T_return

    def _clamp(self, T):
        return min(max(T, self.fluid.T_min), self.fluid.T_max)

    # -- closure ---------------------------------------------------------------
    def _bracket(self, n=200):
        lo, hi = self.fluid.T_min, self.fluid.T_max
        vals = [(lo + (hi - lo) * i / n, self.residual(lo + (hi - lo) * i / n))
                for i in range(n + 1)]
        for (Ta, ra), (Tb, rb) in zip(vals, vals[1:]):
            if ra == 0.0:
                return Ta, Ta
            if ra * rb < 0.0:
                return Ta, Tb
        # The residual decreases with T_return, so a negative value at the freezing point
        # means the root lies below it, and a positive value at the boiling point above.
        direction = (f"below the {lo:.1f} K freezing limit" if vals[0][1] < 0.0
                     else f"above the {hi:.1f} K boiling limit")
        raise LoopDomainError(
            f"no FLiBe loop solution inside {lo:.1f}-{hi:.1f} K: the return temperature "
            f"would lie {direction}. Reduce UA or raise the secondary inlet temperature."
        )

    def solve(self, max_iter=200, tol=1e-9):
        a, b = self._bracket()
        T_return = (a if a == b else
                    brentq(self.residual, a, b, maxiter=max_iter, xtol=1e-12, rtol=8.9e-16))
        r = self.residual(T_return)
        if abs(r) > tol:
            raise RuntimeError(f"loop closure did not converge: residual {r:.3e} K > {tol:.1e}")

        hx = self._exchanger(T_return)
        cold = self.fluid.state_TP(T_return, self.P)
        hot = self._core_outlet(T_return)
        W_pump = self.pump_work(T_return)
        return dict(
            T_return=T_return,
            T_core_out=hot.T,
            T_rise=hot.T - T_return,
            T_sec_out=hx.cold_out.T,
            Q_core=self.Q_core,
            Q_hx=hx.q,
            W_pump=W_pump,
            eps=hx.eps,
            NTU=hx.NTU,
            Cr=hx.Cr,
            C_hot=hx.C_hot,
            C_cold=hx.C_cold,
            residual=r,
            balance=abs(self.Q_core + W_pump - hx.q) / self.Q_core,
            freeze_margin=T_return - self.fluid.T_min,
            states=(cold, hot, hx.hot_out, hx.cold_out),
        )
