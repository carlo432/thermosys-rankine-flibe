"""Coupled FLiBe -> reheat steam plant."""
import warnings
from unittest import mock

import pytest

from thermosys.cycles import CoupledPlant, LoopDomainError, ReheatCycle
from thermosys.fluids import FlibeFluid


def make_plant(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return CoupledPlant(FlibeFluid(), 3900.0, 500e6, **kw)


def solve(**kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_plant(**kw).solve()


def test_global_energy_balance():
    """Q_to_steam = Q_core + W_pump: pump work is added to the coolant, not removed."""
    r = solve()
    assert r["balance"] < 1e-9
    assert r["Q_to_steam"] == pytest.approx(r["Q_core"] + r["W_pump_flibe"], rel=1e-12)


def test_closure_converged():
    assert abs(solve()["residual"]) < 1e-9


@pytest.mark.parametrize("max_iter", [20, 200, 2000])
def test_invariant_to_iteration_cap(max_iter):
    base = solve()["T_return"]
    assert make_plant().solve(max_iter=max_iter)["T_return"] == pytest.approx(base, rel=1e-10)


def test_convergence_guard_rejects_a_bad_root():
    with mock.patch("thermosys.cycles.plant.brentq", return_value=1500.0):
        with pytest.raises(RuntimeError, match="did not converge"):
            make_plant().solve()


def test_no_temperature_crossing_at_the_interface():
    r = solve()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        steam = ReheatCycle().solve()
    assert r["approach_K"] > 0.0
    assert r["T_core_out"] > steam["states"][2].T


def test_plant_efficiency_is_below_the_standalone_cycle():
    r = solve()
    assert r["eta_plant"] < r["eta_steam_alone"]
    assert r["penalty_points"] > 0.0


def test_flibe_stays_liquid():
    r = solve()
    f = FlibeFluid()
    assert f.T_min < r["T_return"] < f.T_max
    assert f.T_min < r["T_core_out"] < f.T_max
    assert r["freeze_margin"] > 0.0


def test_steam_flow_absorbs_the_duty():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        steam = ReheatCycle().solve()
    r = solve()
    # q_in, not q_boiler: the primary supplies the reheater too (see CoupledPlant docstring)
    assert r["mdot_steam"] * steam["q_in"] == pytest.approx(r["Q_to_steam"], rel=1e-9)


@pytest.mark.parametrize("UA,T_expected", [(2.0e6, 820.6190), (3.0e6, 759.3162)])
def test_return_temperature_responds_to_UA(UA, T_expected):
    """R5 again: bypassing the exchanger cannot reproduce this."""
    assert solve(UA=UA)["T_return"] == pytest.approx(T_expected, abs=1e-3)
    assert solve(UA=UA)["T_return"] != pytest.approx(solve()["T_return"], abs=1e-3)


@pytest.mark.parametrize("UA", [1.0e4, 1.2e7])
def test_infeasible_coupling_raises(UA):
    """Too little UA and the primary cannot reject its heat below the boiling point; too
    much and the return temperature falls through the freezing point."""
    with pytest.raises(LoopDomainError):
        solve(UA=UA)
