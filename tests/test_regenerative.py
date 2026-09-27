"""Regenerative cycle with one open FWH."""
import warnings

import pytest

from thermosys.cycles import RegenerativeCycle, RankineCycle


def solve(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return RegenerativeCycle(**kw).solve()


def test_ideal_regenerative():
    r = solve()
    assert r["y"] == pytest.approx(0.214965, abs=1e-6)
    assert r["eta_th"] == pytest.approx(0.421198, abs=1e-6)
    assert r["T_fwh_sat"] == pytest.approx(453.0280, abs=1e-4)
    assert r["q_in"] == pytest.approx(2629101.8, rel=1e-6)
    assert r["w_turbine"] == pytest.approx(1116032.3, rel=1e-6)
    assert r["w_pump"] == pytest.approx(8659.7481, rel=1e-6)
    for key, want in (("h1", 191805.9), ("h2", 192805.9), ("h3", 762515.1),
                      ("h4", 770389.8), ("h5", 3399491.7), ("h6", 2843047.6),
                      ("h7", 2130228.1)):
        assert r[key] == pytest.approx(want, rel=1e-6)


def test_extraction_fraction_is_solved_not_assumed():
    """R2: y satisfies the FWH energy balance it was derived from."""
    assert solve()["fwh_residual"] < 1e-12


def test_fwh_outlet_is_saturated_liquid():
    r = solve()
    s3 = r["states"][2]
    assert s3.x == pytest.approx(0.0, abs=1e-12)
    assert s3.P == pytest.approx(1e6)


def test_regeneration_beats_the_plain_cycle():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plain = RankineCycle(P_boiler=8e6, P_cond=1e4, T_boiler=773.15).solve()
    assert solve()["eta_th"] > plain["eta_th"]


def test_energy_closure():
    r = solve()
    assert abs(r["q_in"] - r["q_out"] - r["w_net"]) / r["q_in"] < 1e-9


@pytest.mark.parametrize("P_extract", [0.3e6, 2e6])
def test_efficiency_responds_to_extraction_pressure(P_extract):
    assert solve(P_extract=P_extract)["eta_th"] != pytest.approx(solve()["eta_th"], abs=1e-5)
