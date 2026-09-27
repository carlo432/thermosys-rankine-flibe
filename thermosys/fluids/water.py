"""CoolProp-backed water/steam.

This is the ONLY module in the package permitted to import CoolProp.
Components reach properties through the Fluid protocol so that FLiBe is a swap.
"""
import CoolProp.CoolProp as CP

from thermosys.state import State

_FLUID = "Water"


def _quality(P, h):
    """CoolProp reports Q = -1 outside the dome; the State contract wants None there."""
    try:
        x = CP.PropsSI("Q", "P", P, "H", h, _FLUID)
    except Exception:
        return None
    return None if x < 0.0 or x > 1.0 else x


class WaterFluid:
    name = "water"
    T_min = 273.16
    T_max = 2000.0

    def _state(self, T, P, h, s):
        return State(self, T, P, h, s, _quality(P, h))

    def state_TP(self, T, P):
        return self._state(T, P,
                           CP.PropsSI("H", "T", T, "P", P, _FLUID),
                           CP.PropsSI("S", "T", T, "P", P, _FLUID))

    def state_PH(self, P, h):
        return self._state(CP.PropsSI("T", "P", P, "H", h, _FLUID), P, h,
                           CP.PropsSI("S", "P", P, "H", h, _FLUID))

    def state_PS(self, P, s):
        return self._state(CP.PropsSI("T", "P", P, "S", s, _FLUID), P,
                           CP.PropsSI("H", "P", P, "S", s, _FLUID), s)

    def state_PQ(self, P, Q):
        h = CP.PropsSI("H", "P", P, "Q", Q, _FLUID)
        return State(self, CP.PropsSI("T", "P", P, "Q", Q, _FLUID), P, h,
                     CP.PropsSI("S", "P", P, "Q", Q, _FLUID), Q)

    def rho(self, T, P):
        return CP.PropsSI("D", "T", T, "P", P, _FLUID)

    def cp(self, T, P):
        return CP.PropsSI("C", "T", T, "P", P, _FLUID)

    def mu(self, T, P):
        return CP.PropsSI("V", "T", T, "P", P, _FLUID)

    def k(self, T, P):
        return CP.PropsSI("L", "T", T, "P", P, _FLUID)

    def Pr(self, T, P):
        return CP.PropsSI("Prandtl", "T", T, "P", P, _FLUID)
