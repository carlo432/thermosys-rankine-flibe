"""FLiBe loop closure. Falsifiability rules R2-R5 are all exercised here."""
import math
import warnings
from unittest import mock

import pytest

from thermosys.components.hx import counterflow_hx, effectiveness_counterflow
from thermosys.cycles import FlibeLoop, LoopDomainError
from thermosys.fluids import ConstantCpFluid, FlibeFluid


def make_loop(UA=4.5e6, T_sec_in=800.0, **kw):
    return FlibeLoop(FlibeFluid(), 3900.0, 500e6, ConstantCpFluid(cp=5200.0),
                     2500.0, T_sec_in, UA, **kw)


def test_loop_anchor():
    r = make_loop().solve()
    assert r["C_hot"] == pytest.approx(9305400.0, rel=1e-12)
    assert r["C_cold"] == pytest.approx(13000000.0, rel=1e-12)
    assert r["Cr"] == pytest.approx(0.715800000, rel=1e-9)
    assert r["NTU"] == pytest.approx(0.483590173, rel=1e-8)
    assert r["eps"] == pytest.approx(0.341411051, rel=1e-8)
    assert r["T_return"] == pytest.approx(903.691595, abs=1e-5)
    assert r["T_core_out"] == pytest.approx(957.445088, abs=1e-5)
    assert r["T_rise"] == pytest.approx(53.753493, abs=1e-5)
    assert r["W_pump"] == pytest.approx(197754.677125, rel=1e-9)
    assert r["Q_hx"] == pytest.approx(500197754.677125, rel=1e-12)
    assert r["T_sec_out"] == pytest.approx(838.476750, abs=1e-5)
    assert r["freeze_margin"] == pytest.approx(171.691595, abs=1e-5)


def test_energy_balance_is_a_real_check():
    """Q_hx comes from the exchanger, not from Q_core - W_pump (rule R2)."""
    r = make_loop().solve()
    assert r["balance"] < 1e-9


@pytest.mark.parametrize("UA,T_expected", [(3.0e6, 959.211722),
                                           (4.5e6, 903.691595),
                                           (6.0e6, 875.960871)])
def test_UA_sensitivity(UA, T_expected):
    """R5: a build that bypasses the exchanger cannot reproduce this table."""
    assert make_loop(UA=UA).solve()["T_return"] == pytest.approx(T_expected, abs=1e-5)


@pytest.mark.parametrize("max_iter", [20, 200, 2000])
def test_answer_is_invariant_to_iteration_cap(max_iter):
    """R4: an unconverged fixed-point iteration must not be reported as a solution."""
    assert make_loop().solve(max_iter=max_iter)["T_return"] == pytest.approx(
        903.691595120, rel=1e-10)


def test_convergence_guard_rejects_a_bad_root():
    """solve() must verify the root rather than trust the solver: hand it a wrong answer
    and it has to notice. Without this the guard is unfalsified defensive code."""
    with mock.patch("thermosys.cycles.flibe_loop.brentq", return_value=1200.0):
        with pytest.raises(RuntimeError, match="did not converge"):
            make_loop().solve()


def test_out_of_domain_raises_named_error():
    """§5.2: cold secondary + huge UA drives the root below the freezing point."""
    with pytest.raises(LoopDomainError, match="freezing"):
        make_loop(UA=4.0e7, T_sec_in=600.0).solve()


@pytest.mark.parametrize("NTU,Cr,eps", [
    (0.5, 0.00, 0.393469340), (0.5, 0.50, 0.362265573),
    (1.0, 0.00, 0.632120559), (1.0, 0.50, 0.564733402), (1.0, 1.00, 0.500000000),
    (2.0, 0.50, 0.774600326), (2.0, 1.00, 0.666666667),
    (3.0, 0.25, 0.918811274), (5.0, 0.80, 0.895740033), (10.0, 1.00, 0.909090909),
])
def test_effectiveness_table(NTU, Cr, eps):
    assert effectiveness_counterflow(NTU, Cr) == pytest.approx(eps, abs=1e-9)


def test_effectiveness_limits():
    assert effectiveness_counterflow(1.0, 0.0) == pytest.approx(1 - math.exp(-1), abs=1e-12)
    assert effectiveness_counterflow(10.0, 1.0) == pytest.approx(10 / 11, abs=1e-12)
    assert effectiveness_counterflow(0.0, 0.5) == pytest.approx(0.0, abs=1e-15)


def test_no_temperature_crossing():
    f, sec = FlibeFluid(), ConstantCpFluid(cp=5200.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cold_in = sec.state_TP(800.0, 1e5)
        hx = counterflow_hx(f.state_TP(1000.0, 1e5), cold_in, 3900.0, 2500.0, 4.5e6)
    assert hx.hot_out.T >= cold_in.T
    assert hx.cold_out.T <= 1000.0
    assert hx.q > 0.0
