"""Generic constant-property secondary stream.

Not a third working fluid: the FLiBe loop specification rejects heat to a
secondary stream characterised only by (mdot, cp, T_in), and `counterflow_hx` takes
States.  This is the minimal Fluid that lets that stream be expressed in the existing
component signatures.
"""
import math

from thermosys.state import State


class ConstantCpFluid:
    """Incompressible, constant-cp stream. h = cp*(T - T_ref), s = cp*ln(T/T_ref)."""

    def __init__(self, cp=5200.0, rho=1000.0, mu=1.0e-3, k=0.6,
                 T_ref=273.15, T_min=200.0, T_max=3000.0, name="secondary"):
        self.name = name
        self._cp, self._rho, self._mu, self._k = cp, rho, mu, k
        self.T_ref, self.T_min, self.T_max = T_ref, T_min, T_max

    def state_TP(self, T, P):
        return State(self, T, P, self._cp * (T - self.T_ref),
                     self._cp * math.log(T / self.T_ref), None)

    def state_PH(self, P, h):
        return self.state_TP(self.T_ref + h / self._cp, P)

    def state_PS(self, P, s):
        return self.state_TP(self.T_ref * math.exp(s / self._cp), P)

    def state_PQ(self, P, Q):
        raise NotImplementedError("constant-cp stream has no two-phase quality")

    def rho(self, T, P=0.0):
        return self._rho

    def cp(self, T, P=0.0):
        return self._cp

    def mu(self, T, P=0.0):
        return self._mu

    def k(self, T, P=0.0):
        return self._k

    def Pr(self, T, P=0.0):
        return self._mu * self._cp / self._k
