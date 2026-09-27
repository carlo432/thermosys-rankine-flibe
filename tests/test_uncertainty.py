"""FLiBe property uncertainty propagation, and the correlations that make it non-vacuous."""
import warnings

import pytest

from thermosys.analysis.uncertainty import (DEFAULT_BANDS, corner_study, one_at_a_time,
                                            solve_with_computed_UA)
from thermosys.components.htc import dittus_boelter, reynolds
from thermosys.fluids import FlibeFluid
from thermosys.fluids.perturbed import PerturbedFluid


def test_perturbed_fluid_scales_only_what_it_is_told_to():
    base = FlibeFluid()
    f = PerturbedFluid(base, mu=1.2, k=0.85)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert f.mu(900.0) == pytest.approx(1.2 * base.mu(900.0), rel=1e-12)
        assert f.k(900.0) == pytest.approx(0.85 * base.k(900.0), rel=1e-12)
        assert f.rho(900.0) == pytest.approx(base.rho(900.0), rel=1e-12)
        assert f.cp(900.0) == pytest.approx(base.cp(900.0), rel=1e-12)
        # Pr must follow the perturbed properties, not the base ones
        assert f.Pr(900.0) == pytest.approx(1.2 / 0.85 * base.Pr(900.0), rel=1e-12)


def test_dittus_boelter_is_turbulent_and_in_range():
    f = FlibeFluid()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Re = reynolds(f, 900.0, 1e5, 3900.0 / 2000, 0.02)
        Pr = f.Pr(900.0, 1e5)
        h = dittus_boelter(f, 900.0, 1e5, 3900.0 / 2000, 0.02)
    assert Re > 10000.0, "correlation requires turbulent flow"
    assert 0.6 < Pr < 160.0
    assert h > 0.0


def test_UA_is_computed_from_properties_not_specified():
    """If UA did not depend on mu and k, the uncertainty study below would be vacuous."""
    base = FlibeFluid()
    nominal = solve_with_computed_UA(PerturbedFluid(base))["UA"]
    assert solve_with_computed_UA(PerturbedFluid(base, mu=1.2))["UA"] != pytest.approx(
        nominal, rel=1e-6)
    assert solve_with_computed_UA(PerturbedFluid(base, k=1.15))["UA"] != pytest.approx(
        nominal, rel=1e-6)


def test_conductivity_and_viscosity_dominate():
    """The ranking is the point: the two correlations extrapolated beyond their fitted
    ranges are the two that move the answer."""
    o = one_at_a_time(FlibeFluid())
    e = o["effects"]
    assert e["k"]["range_K"] > 10.0
    assert e["mu"]["range_K"] > 10.0
    assert e["cp"]["range_K"] < 1.0
    assert e["rho"]["range_K"] < 1.0


def test_every_corner_solves_and_stays_liquid():
    rows = corner_study(FlibeFluid())
    solved = [r for r in rows if "T_return" in r]
    assert len(solved) == len(rows), [r for r in rows if "error" in r]
    assert min(r["freeze_margin"] for r in solved) > 50.0


def test_total_spread_is_reported_and_bounded():
    rows = [r for r in corner_study(FlibeFluid()) if "T_return" in r]
    Ts = [r["T_return"] for r in rows]
    spread = max(Ts) - min(Ts)
    assert 20.0 < spread < 60.0, f"property uncertainty moves T_return by {spread:.1f} K"


def test_bands_are_the_documented_ones():
    assert DEFAULT_BANDS == dict(mu=0.20, k=0.15, rho=0.02, cp=0.03)
