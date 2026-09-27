"""Simple Rankine cases and the cycle invariants."""
import warnings

import pytest

from thermosys.cycles import RankineCycle

CASE_A = dict()
CASE_B = dict(P_boiler=8e6, P_cond=1e4, T_boiler=773.15, eta_p=0.85, eta_t=0.88)
CASE_C = dict(P_boiler=8e6, P_cond=1e4, T_boiler=773.15)


def solve(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return RankineCycle(**kw).solve()


@pytest.mark.parametrize("kw,eta,x4", [(CASE_A, 0.248878, 0.796538),
                                       (CASE_B, 0.346279, 0.874033),
                                       (CASE_C, 0.394173, 0.810359)])
def test_case_efficiency_and_quality(kw, eta, x4):
    r = solve(**kw)
    assert r["eta_th"] == pytest.approx(eta, abs=1e-6)
    assert r["states"][3].x == pytest.approx(x4, abs=1e-6)


def test_case_a_enthalpies():
    r = solve(**CASE_A)
    for key, want in (("h1", 384437.0), ("h2", 387469.0), ("h3", 2803150.0),
                      ("h4", 2198910.0)):
        assert r[key] == pytest.approx(want, rel=1e-5)
    assert r["q_in"] == pytest.approx(2415680.0, rel=1e-5)
    assert r["w_turbine"] == pytest.approx(604244.0, rel=1e-5)
    assert r["w_pump"] == pytest.approx(3031.98, rel=1e-5)
    assert r["back_work_ratio"] == pytest.approx(0.00501781, rel=1e-5)


@pytest.mark.parametrize("kw", [CASE_A, CASE_B, CASE_C])
def test_energy_closure(kw):
    r = solve(**kw)
    assert abs(r["q_in"] - r["q_out"] - r["w_net"]) / r["q_in"] < 1e-9


@pytest.mark.parametrize("kw", [CASE_A, CASE_C])
def test_isentropic_recovery(kw):
    """eta_s = 1 must reproduce the isentropic outlet. The pump path floors at 2.7e-9
    relative (CoolProp state_PS -> state_PH round trip), so 1e-8 rather than the originally targeted
    1e-9; the turbine path is exact."""
    s1, s2, s3, s4 = solve(**kw)["states"]
    assert s2.s == pytest.approx(s1.s, rel=1e-8)
    assert s4.s == pytest.approx(s3.s, rel=1e-12)


def test_entropy_does_not_decrease():
    s1, s2, s3, s4 = solve(**CASE_B)["states"]
    assert s2.s > s1.s
    assert s4.s > s3.s


def test_superheat_and_pressure_raise_efficiency():
    assert solve(**CASE_C)["eta_th"] > solve(**CASE_A)["eta_th"]


def test_pump_work_matches_incompressible_estimate():
    """Catches the classic isentropic-direction error: v*dP/eta vs the enthalpy path."""
    r = solve(**CASE_B)
    s1 = r["states"][0]
    # evaluated 0.1 K subcooled: a (T, P) property call exactly on the saturation line is
    # ambiguous and CoolProp may return the vapour branch
    v1 = 1.0 / s1.fluid.rho(s1.T - 0.1, s1.P)
    assert r["w_pump"] == pytest.approx(9480.47, rel=5e-3)
    assert r["w_pump"] == pytest.approx(v1 * (8e6 - 1e4) / 0.85, rel=5e-3)


def test_wetness_warning_fires_only_when_wet():
    with pytest.warns(UserWarning, match="turbine exit quality"):
        RankineCycle(**CASE_A).solve()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        RankineCycle(P_boiler=1e5, P_cond=1e4, T_boiler=773.15).solve()
