"""Moving-boundary boiler."""
import numpy as np
import pytest

from thermosys.fluids import WaterFluid
from thermosys.transient.boiler_zones import MovingBoundaryBoiler

P = 8e6


def make(q_total=1.4e8, n_cells=60, **kw):
    f = WaterFluid()
    h_in = f.state_PQ(P, 0.0).h - 3.0e5          # 300 kJ/kg subcooled
    return MovingBoundaryBoiler(fluid=f, P=P, mdot=50.0, h_in=h_in,
                                q_total=q_total, length=20.0, M_total=500.0,
                                n_cells=n_cells, **kw)


def test_saturation_enthalpies_come_from_the_fluid():
    b = make()
    assert b.h_f == pytest.approx(1317311.05, rel=1e-6)
    assert b.h_g == pytest.approx(2758683.20, rel=1e-6)


@pytest.mark.parametrize("n_cells", [10, 20, 60, 200])
def test_zone_diagnosis_is_grid_independent(n_cells):
    """Interpolating on the downstream-face grid makes the diagnosis exact for the linear
    steady profile, so refining the mesh must not move the boundaries at all."""
    b = make(n_cells=n_cells)
    assert b.zone_lengths(b.steady_profile()) == pytest.approx(
        b.analytic_zone_lengths(), abs=1e-9)


def test_analytic_zone_lengths():
    b = make()
    L1, L2, L3 = b.analytic_zone_lengths()
    q_prime = b.q_total / b.length
    assert L1 == pytest.approx(b.mdot * (b.h_f - b.h_in) / q_prime, rel=1e-12)
    assert L2 == pytest.approx(b.mdot * (b.h_g - b.h_f) / q_prime, rel=1e-12)
    assert L1 + L2 + L3 == pytest.approx(b.length, rel=1e-12)


def test_analytic_zone_lengths_are_clipped_when_a_zone_overruns_the_tube():
    """At low heat input the unclipped two-phase length is far longer than the tube; the
    closed form has to clip, and must still agree with the diagnosed boundaries."""
    b = make(q_total=2.0e7, n_cells=200)
    q_prime = b.q_total / b.length
    assert b.mdot * (b.h_g - b.h_f) / q_prime > b.length      # unclipped overruns
    L1, L2, L3 = b.analytic_zone_lengths()
    assert (L1, L2, L3) == pytest.approx((15.0, 5.0, 0.0), abs=1e-9)
    assert L1 + L2 + L3 == pytest.approx(b.length, rel=1e-12)
    assert b.zone_lengths(b.steady_profile()) == pytest.approx((L1, L2, L3), abs=1e-9)


def test_dynamic_solution_recovers_the_steady_zones():
    b = make()
    sol = b.solve(t_final=600.0)
    assert b.zone_lengths(sol.y[:, -1]) == pytest.approx(b.analytic_zone_lengths(), abs=1e-6)


def test_outlet_enthalpy_closes_the_energy_balance():
    b = make()
    sol = b.solve(t_final=600.0)
    assert float(sol.h_out[-1]) == pytest.approx(b.h_in + b.q_total / b.mdot, rel=1e-9)


def test_superheat_zone_vanishes_gracefully():
    """Where the classical three-state formulation goes singular, this must just report
    zero length and a two-phase outlet."""
    b = make(q_total=2.0e7, n_cells=200)
    L1, L2, L3 = b.zone_lengths(b.steady_profile())
    assert L3 == pytest.approx(0.0, abs=1e-9)
    assert L1 + L2 == pytest.approx(b.length, rel=1e-12)
    assert 0.0 < b.outlet_quality(b.steady_profile()) < 1.0


def test_dryout_moves_upstream_as_power_rises():
    lengths = []
    for q in (1.0e8, 1.4e8, 2.0e8):
        b = make(q_total=q, n_cells=200)
        L1, L2, _ = b.zone_lengths(b.steady_profile())
        lengths.append(L1 + L2)
    assert lengths[0] > lengths[1] > lengths[2]


def test_transport_delay_matches_tube_fill_time():
    """A distributed model must stay silent for the fill time; a lumped CV cannot."""
    b = make(n_cells=200)
    steady = b.steady_profile()
    stepped = make(n_cells=200)
    stepped.h_in = b.h_in + 2.0e5
    sol = stepped.solve(h0=steady, t_final=60.0)
    tau = b.M_total / b.mdot                     # 10 s
    frac = lambda t: (float(sol.sol(t)[-1]) - steady[-1]) / 2.0e5
    assert frac(0.5 * tau) < 0.02                # quiet well inside the fill time
    assert 0.3 < frac(tau) < 0.7                 # arriving at the fill time
    assert frac(3 * tau) > 0.98                  # complete afterwards


def test_zone_lengths_never_exceed_the_tube():
    for q in (1.0e6, 2.0e7, 1.4e8, 1.0e9):
        b = make(q_total=q, n_cells=100)
        L = b.zone_lengths(b.steady_profile())
        assert all(x >= -1e-12 for x in L)
        assert sum(L) == pytest.approx(b.length, rel=1e-12)
