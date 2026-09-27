"""Reheat cycle."""
import warnings

import pytest

from thermosys.cycles import ReheatCycle, RankineCycle


def solve(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ReheatCycle(**kw).solve()


def test_ideal_reheat():
    r = solve()
    assert r["eta_th"] == pytest.approx(0.408355, abs=1e-6)
    assert r["x6"] == pytest.approx(0.948708, abs=1e-6)
    assert r["q_in"] == pytest.approx(3835680.1, rel=1e-6)
    assert r["w_turbine"] == pytest.approx(1574379.4, rel=1e-6)
    assert r["w_pump"] == pytest.approx(8058.399, rel=1e-6)
    for key, want in (("h1", 191805.9), ("h2", 199864.3), ("h3", 3399491.7),
                      ("h4", 2843047.6), ("h5", 3479100.4), ("h6", 2461165.1)):
        assert r[key] == pytest.approx(want, rel=1e-6)


def test_non_ideal_reheat():
    r = solve(eta_p=0.85, eta_t=0.90)
    assert r["eta_th"] == pytest.approx(0.372481, abs=1e-6)
    assert r["x6"] == pytest.approx(0.991263, abs=1e-6)
    assert r["q_in"] == pytest.approx(3778613.7, rel=1e-6)
    assert r["w_turbine"] == pytest.approx(1416941.5, rel=1e-6)
    assert r["w_pump"] == pytest.approx(9480.470, rel=1e-6)


def test_reheat_beats_the_plain_cycle():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plain = RankineCycle(P_boiler=8e6, P_cond=1e4, T_boiler=773.15).solve()
    assert solve()["eta_th"] > plain["eta_th"]


def test_reheat_cures_the_wetness_problem():
    """The engineering point of reheat, as a test rather than a comment."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plain = RankineCycle(P_boiler=8e6, P_cond=1e4, T_boiler=773.15).solve()
    assert plain["states"][3].x < 0.88
    assert solve()["x6"] > 0.88
    with warnings.catch_warnings():          # and therefore fires no wetness warning
        warnings.simplefilter("error")
        ReheatCycle().solve()


def test_energy_closure():
    for kw in ({}, dict(eta_p=0.85, eta_t=0.90)):
        r = solve(**kw)
        assert abs(r["q_in"] - r["q_out"] - r["w_net"]) / r["q_in"] < 1e-9


@pytest.mark.parametrize("P_reheat", [0.5e6, 2e6, 4e6])
def test_efficiency_responds_to_reheat_pressure(P_reheat):
    """R5: the answer must depend on the governing input."""
    assert solve(P_reheat=P_reheat)["eta_th"] != pytest.approx(solve()["eta_th"], abs=1e-5)


def test_reheat_pressure_optimum_is_interior():
    etas = {P: solve(P_reheat=P)["eta_th"] for P in (0.5e6, 1e6, 2e6, 4e6)}
    best = max(etas, key=etas.get)
    assert best not in (0.5e6, 4e6), f"expected an interior optimum, got {etas}"
