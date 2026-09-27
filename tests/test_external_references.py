"""Checks against third-party published solutions, not against my own anchors.

Every other test in this suite compares the code to reference values I generated myself.
These four compare it to numbers published by someone else, which is a different and
stronger claim. Tolerances are set by the publications' own rounding.
"""
import warnings

import pytest

from thermosys.cycles import RankineCycle, ReheatCycle


def solve_reheat(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ReheatCycle(**kw).solve()


def solve_plain(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return RankineCycle(**kw).solve()


def test_published_reheat_enthalpies_8MPa_4MPa_10kPa():
    """Published worked solution: ideal reheat between 8 MPa and 10 kPa, reheat at 4 MPa,
    500 C at both turbine inlets -> HP exit 3185 kJ/kg, LP exit 2247 kJ/kg (4 s.f.)."""
    r = solve_reheat(P_boiler=8e6, P_reheat=4e6, P_cond=1e4,
                     T_boiler=773.15, T_reheat=773.15)
    assert r["h4"] / 1000.0 == pytest.approx(3185.0, abs=0.5)
    assert r["h6"] / 1000.0 == pytest.approx(2247.0, abs=0.5)


def test_published_quality_15MPa_600C_plain():
    """Ohio University applied thermo notes, ideal Rankine 15 MPa / 600 C / 10 kPa:
    LP turbine exit quality quoted as 80%."""
    r = solve_plain(P_boiler=15e6, P_cond=1e4, T_boiler=873.15)
    assert r["states"][3].x == pytest.approx(0.80, abs=0.005)


def test_published_quality_15MPa_600C_reheat():
    """Same source, with reheat from 15 MPa to 1 MPa and back to 600 C: quality 98%."""
    r = solve_reheat(P_boiler=15e6, P_reheat=1e6, P_cond=1e4,
                     T_boiler=873.15, T_reheat=873.15)
    assert r["x6"] == pytest.approx(0.98, abs=0.005)


def test_reheat_raises_quality_at_the_published_conditions():
    plain = solve_plain(P_boiler=15e6, P_cond=1e4, T_boiler=873.15)
    reheat = solve_reheat(P_boiler=15e6, P_reheat=1e6, P_cond=1e4,
                          T_boiler=873.15, T_reheat=873.15)
    assert plain["states"][3].x < 0.85 < reheat["x6"]
    assert reheat["eta_th"] > plain["eta_th"]
