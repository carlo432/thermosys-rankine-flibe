"""Convective heat-transfer correlations.

These exist so that UA can be *computed from properties* rather than specified. Without
them the loop's answer is insensitive to viscosity and conductivity, and any claim about
FLiBe property uncertainty would be vacuous -- mu and k would not enter the model at all.
"""
import math


def reynolds(fluid, T, P, mdot, D):
    """Re = 4*mdot / (pi*D*mu) for a round tube."""
    return 4.0 * mdot / (math.pi * D * fluid.mu(T, P))


def dittus_boelter(fluid, T, P, mdot, D, heating=True):
    """h = 0.023 Re^0.8 Pr^n k/D, n = 0.4 heating / 0.3 cooling.

    Validity, stated because FLiBe sits near the edge of it: Re > 10000, 0.6 < Pr < 160,
    L/D > 10. FLiBe at 900 K has Pr ~ 14, comfortably inside; a cold-leg excursion toward
    the freezing point pushes Pr past 35 and the correlation degrades.
    """
    Re = reynolds(fluid, T, P, mdot, D)
    Pr = fluid.Pr(T, P)
    Nu = 0.023 * Re ** 0.8 * Pr ** (0.4 if heating else 0.3)
    return Nu * fluid.k(T, P) / D


def ua_from_htc(h, area):
    return h * area


def loop_UA(fluid, T, P, mdot, D, area, n_tubes=1, heating=False):
    """Overall UA taking the salt-side film as controlling (a documented simplification:
    the secondary film and wall resistance are neglected, so this is an upper bound)."""
    h = dittus_boelter(fluid, T, P, mdot / n_tubes, D, heating=heating)
    return ua_from_htc(h, area)
