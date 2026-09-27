"""Water properties and the single-phase quality contract."""
import pytest

from thermosys.fluids import WaterFluid

SAT = [
    (1e4, 318.9563, 191805.94, 2583858.67),
    (75e3, 364.9080, 384437.44, 2662386.22),
    (101325.0, 373.1243, 419057.73, 2675529.33),
    (3e6, 507.0031, 1008344.61, 2803153.10),
    (8e6, 568.1577, 1317311.05, 2758683.20),
]


@pytest.mark.parametrize("P,T_sat,h_f,h_g", SAT)
def test_saturation_table(P, T_sat, h_f, h_g):
    f = WaterFluid()
    assert f.state_PQ(P, 0).T == pytest.approx(T_sat, rel=1e-6)
    assert f.state_PQ(P, 0).h == pytest.approx(h_f, rel=1e-6)
    assert f.state_PQ(P, 1).h == pytest.approx(h_g, rel=1e-6)


def test_quality_is_none_outside_the_dome():
    """CoolProp reports Q = -1 in single phase; the State contract says None."""
    f = WaterFluid()
    assert f.state_TP(300.0, 8e6).x is None       # compressed liquid
    assert f.state_TP(700.0, 1e3).x is None       # superheated vapour
    assert f.state_PQ(1e4, 0.5).x == pytest.approx(0.5)


def test_round_trip_consistency():
    f = WaterFluid()
    a = f.state_TP(500.0, 3e6)
    assert f.state_PH(3e6, a.h).T == pytest.approx(500.0, rel=1e-9)
    assert f.state_PS(3e6, a.s).h == pytest.approx(a.h, rel=1e-9)
