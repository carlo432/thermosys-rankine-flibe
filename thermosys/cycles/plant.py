"""Coupled plant: FLiBe primary loop driving the reheat steam cycle.

The FLiBe loop's secondary stream IS the steam cycle boiler.  No new components: the
interface is the same counterflow exchanger, with the steam side represented by its
effective capacity rate

    C_c = mdot_steam * (h3 - h2) / (T3 - T2)

which is the standard mean-cp linearisation of a boiling stream.  It is a linearisation,
not an exact two-phase HX model, and is documented as such.

One unknown (the FLiBe return temperature) closes everything: the pump work follows from
it, the steam flow follows from the duty, and the exchanger then has to reproduce the
return temperature it started from.
"""
from scipy.optimize import brentq

from thermosys.components import heat_source
from thermosys.components.hx import effectiveness_counterflow
from thermosys.cycles.flibe_loop import LoopDomainError
from thermosys.cycles.reheat import ReheatCycle


class CoupledPlant:
    def __init__(self, flibe, mdot_flibe, Q_core, steam_cycle=None, UA=2.5e6,
                 deltaP=1.0e5, P=1.0e5):
        self.flibe = flibe
        self.mdot = mdot_flibe
        self.Q_core = Q_core
        self.steam = steam_cycle or ReheatCycle()
        self.UA = UA
        self.deltaP = deltaP
        self.P = P
        self._steam_result = self.steam.solve()

    # -- steam side, fixed by its own boundary conditions ----------------------
    @property
    def _steam_interface(self):
        """Feedwater state, live-steam state, and the TOTAL specific heat the primary has
        to supply.  It is q_in, not q_boiler: in a FLiBe-heated plant the reheater is fed
        from the same primary loop, so charging only the boiler duty would hand the cycle
        its reheat energy for free and push plant efficiency above the standalone cycle."""
        r = self._steam_result
        return r["states"][1], r["states"][2], r["q_in"]

    def _duty_and_flow(self, T_return):
        """Pump work, total heat to reject, and the steam flow that absorbs it."""
        W_pump = self.mdot * self.deltaP / self.flibe.rho(self._clamp(T_return), self.P)
        Q_total = self.Q_core + W_pump
        _, _, q_in = self._steam_interface
        return W_pump, Q_total, Q_total / q_in

    def _clamp(self, T):
        return min(max(T, self.flibe.T_min), self.flibe.T_max)

    def residual(self, T_return):
        s2, s3, q_in = self._steam_interface
        W_pump, Q_total, mdot_steam = self._duty_and_flow(T_return)
        C_h = self.mdot * self.flibe.cp(self._clamp(T_return), self.P)
        C_c = mdot_steam * q_in / (s3.T - s2.T)
        T_hot = T_return + Q_total / C_h
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        eps = effectiveness_counterflow(self.UA / C_min, C_min / C_max)
        q = eps * C_min * (T_hot - s2.T)
        return T_hot - q / C_h - T_return

    def _bracket(self, n=300):
        lo, hi = self.flibe.T_min, self.flibe.T_max
        vals = [(lo + (hi - lo) * i / n, self.residual(lo + (hi - lo) * i / n))
                for i in range(n + 1)]
        for (Ta, ra), (Tb, rb) in zip(vals, vals[1:]):
            if ra * rb < 0.0:
                return Ta, Tb
        direction = (f"below the {lo:.1f} K freezing limit" if vals[0][1] < 0.0
                     else f"above the {hi:.1f} K boiling limit")
        raise LoopDomainError(
            f"no coupled-plant solution inside {lo:.1f}-{hi:.1f} K: the FLiBe return "
            f"temperature would lie {direction}. Adjust UA, Q_core, or the steam "
            f"boundary conditions."
        )

    def solve(self, max_iter=200, tol=1e-9):
        a, b = self._bracket()
        T_return = brentq(self.residual, a, b, maxiter=max_iter, xtol=1e-12, rtol=8.9e-16)
        r = self.residual(T_return)
        if abs(r) > tol:
            raise RuntimeError(f"plant closure did not converge: residual {r:.3e} K")

        s2, s3, q_in = self._steam_interface
        W_pump, Q_total, mdot_steam = self._duty_and_flow(T_return)
        cold = self.flibe.state_TP(T_return, self.P)
        hot = heat_source(cold, self.mdot, Q_total)
        if hot.T <= s3.T:
            raise LoopDomainError(
                f"temperature crossing at the interface: FLiBe enters the boiler at "
                f"{hot.T:.2f} K but steam must leave at {s3.T:.2f} K"
            )

        W_steam = mdot_steam * self._steam_result["w_net"]
        W_net = W_steam - W_pump
        return dict(
            T_return=T_return,
            T_core_out=hot.T,
            approach_K=hot.T - s3.T,
            mdot_steam=mdot_steam,
            Q_core=self.Q_core,
            Q_to_steam=Q_total,
            W_pump_flibe=W_pump,
            W_steam=W_steam,
            W_net=W_net,
            eta_plant=W_net / self.Q_core,
            eta_steam_alone=self._steam_result["eta_th"],
            penalty_points=100.0 * (self._steam_result["eta_th"] - W_net / self.Q_core),
            balance=abs(Q_total - (self.Q_core + W_pump)) / self.Q_core,
            residual=r,
            freeze_margin=T_return - self.flibe.T_min,
        )
