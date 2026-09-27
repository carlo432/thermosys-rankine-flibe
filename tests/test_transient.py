"""Transients with real energy balances.

The anchors here are analytic, so deleting any term in a right-hand side breaks a test.
A relaxation model onto a precomputed answer would lack that property.
"""
import warnings

import numpy as np
import pytest

from thermosys.cycles import FlibeLoop
from thermosys.fluids import ConstantCpFluid, FlibeFluid
from thermosys.transient import CVChain, FlibeLoopTransient, SingleCV, erlang_response

CV = dict(M=5000.0, cp=2386.0, mdot=100.0, T_in=800.0, Q=2.0e7)
T_INIT = 750.0


def test_single_cv_constants():
    cv = SingleCV(**CV)
    assert cv.tau == pytest.approx(50.000000000, rel=1e-12)
    assert cv.T_ss == pytest.approx(883.822296731, rel=1e-11)


@pytest.mark.parametrize("t,T", [(0.0, 750.000000000), (10.0, 774.257866950),
                                 (50.0, 834.591824993), (100.0, 865.711418299),
                                 (500.0, 883.816221208)])
def test_single_cv_matches_closed_form(t, T):
    cv = SingleCV(**CV)
    assert cv.analytic(t, T_INIT) == pytest.approx(T, rel=1e-9)
    sol = cv.solve(T_INIT, 500.0, rtol=1e-12, atol=1e-10)
    assert float(sol.sol(t)[0]) == pytest.approx(T, rel=1e-8)


def test_single_cv_responds_to_every_input():
    """R5/R1: each term in the balance changes the answer."""
    base = SingleCV(**CV).T_ss
    assert SingleCV(**{**CV, "Q": 0.0}).T_ss == pytest.approx(CV["T_in"])
    assert SingleCV(**{**CV, "mdot": 200.0}).T_ss < base
    assert SingleCV(**{**CV, "T_in": 900.0}).T_ss > base
    assert SingleCV(**{**CV, "M": 10000.0}).tau == 2 * SingleCV(**CV).tau


# The reference quotes t = 1.0 s to four significant figures (1.721e-04); the exact Erlang
# value is 1.7211562995589347e-04, so that row is checked at rel=1e-3 and the rest at 1e-6.
@pytest.mark.parametrize("t,frac,rel", [(1.0, 1.721e-04, 1e-3), (2.0, 0.003659847, 1e-6),
                                        (5.0, 0.108821981, 1e-6), (10.0, 0.559506715, 1e-6),
                                        (20.0, 0.970747312, 1e-6), (40.0, 0.999983055, 1e-6)])
def test_cv_chain_matches_erlang(t, frac, rel):
    chain = CVChain(N=5, tau_cv=2.0)
    sol = chain.solve(40.0, rtol=1e-12, atol=1e-12)
    assert erlang_response(t, 5, 2.0) == pytest.approx(frac, rel=rel)
    assert float(sol.sol(t)[-1]) == pytest.approx(frac, rel=rel)


def test_chain_produces_transport_delay_a_single_cv_cannot():
    chain = CVChain(N=5, tau_cv=2.0)
    sol = chain.solve(40.0, rtol=1e-12, atol=1e-12)
    assert float(sol.sol(1.0)[-1]) < 1e-3          # still quiet well inside one residence
    single = CVChain(N=1, tau_cv=2.0).solve(40.0, rtol=1e-12, atol=1e-12)
    assert float(single.sol(2.0)[-1]) == pytest.approx(1 - np.exp(-1), rel=1e-6)
    assert float(single.sol(1.0)[-1]) > 100 * float(sol.sol(1.0)[-1])


# --------------------------------------------------------------------------
# loop network
# --------------------------------------------------------------------------
SEC = dict(C_cold=13000000.0, T_sec_in=800.0, UA=4.5e6)


def make_network(**kw):
    base = dict(fluid=FlibeFluid(), mdot=3900.0, Q_core=500e6, **SEC)
    return FlibeLoopTransient(**{**base, **kw})


def make_algebraic(Q=500e6, mdot=3900.0, UA=4.5e6, deltaP=1.0e5):
    return FlibeLoop(FlibeFluid(), mdot, Q, ConstantCpFluid(cp=5200.0), 2500.0, 800.0,
                     UA, deltaP=deltaP)


def test_network_steady_state_matches_algebraic_solve():
    """The check that matters: the ODE network must recover the algebraic solution."""
    net = make_network()
    ref = make_algebraic().solve()
    sol = net.solve(net.steady_guess(870.0), 4000.0)
    assert sol.T_return[-1] == pytest.approx(ref["T_return"], rel=1e-6)
    assert sol.T_core_out[-1] == pytest.approx(ref["T_core_out"], rel=1e-6)


def test_scenario_step_in_core_power():
    """+10% Q must land on the algebraic solution at the new power."""
    net = make_network(Q_core=550e6)
    ref = make_algebraic(Q=550e6).solve()
    sol = net.solve(make_network().steady_guess(903.691595), 4000.0)
    assert sol.T_return[-1] == pytest.approx(ref["T_return"], rel=1e-6)


def test_scenario_step_in_mass_flow():
    """-20% mdot with dP = 0 must raise the core rise by exactly mdot_old/mdot_new."""
    hi = make_network(deltaP=0.0)
    lo = make_network(mdot=3120.0, deltaP=0.0)
    r_hi = hi.solve(hi.steady_guess(900.0), 6000.0)
    r_lo = lo.solve(lo.steady_guess(900.0), 6000.0)
    rise_hi = r_hi.T_core_out[-1] - r_hi.T_return[-1]
    rise_lo = r_lo.T_core_out[-1] - r_lo.T_return[-1]
    assert rise_lo / rise_hi == pytest.approx(3900.0 / 3120.0, rel=1e-6)


def test_scenario_pump_coastdown():
    net = make_network(mdot_func=lambda t: 3900.0 * max(np.exp(-t / 5.0), 0.05))
    sol = net.solve(make_network().steady_guess(903.691595), 20.0, dense_output=True)
    ts = np.linspace(0.0, 20.0, 201)
    core = np.array([sol.sol(t)[0] for t in ts])
    assert np.all(np.diff(core) > -1e-6), "core outlet must rise as flow coasts down"
    peak = float(core.max())
    assert peak > sol.T_core_out[0]
    assert peak < FlibeFluid().T_max


def test_scenario_loss_of_heat_sink_slope_is_analytic():
    """UA -> 0: the loop is adiabatic, so the mean temperature rises at Q/(M_total*cp).
    Impossible to reproduce with a relaxation model."""
    net = make_network(UA_func=lambda t: 0.0, deltaP=0.0)
    y0 = make_network().steady_guess(903.691595)
    sol = net.solve(y0, 100.0)
    expected = net.Q_core / (net.masses.sum() * net.cp)
    slope = (sol.mean_T[-1] - sol.mean_T[0]) / (sol.t[-1] - sol.t[0])
    assert slope == pytest.approx(expected, rel=1e-4)


def test_freeze_excursion_raises():
    """§6.4: loss of heat SOURCE with the sink live drives nodes through 732 K, and that
    is an error, not a comfortable margin."""
    net = make_network(Q_func=lambda t: 0.0, T_sec_in=600.0, UA=6.0e6, deltaP=0.0)
    y0 = np.full(make_network().n_nodes, 760.0)
    with pytest.raises(ValueError, match="froze"):
        net.solve(y0, 3000.0)


def test_integrator_is_stiff_capable():
    from thermosys.transient import integrate
    sol = integrate(lambda t, y: [-1000.0 * y[0]], (0.0, 1.0), [1.0])
    assert sol.success
    assert float(sol.y[0, -1]) == pytest.approx(0.0, abs=1e-6)
