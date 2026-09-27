"""FLiBe property correlations against published reference values (all seven rows, all columns)."""
import warnings

import pytest

from thermosys.fluids import FlibeFluid

# T, rho, mu, nu, Pr(k=1.1 historical), k(T), Pr(k(T))
TABLE = [
    (732.0, 2056.3930, 0.019600, 9.532079e-06, 42.52, 1.2779, 36.60),
    (800.0, 2023.0240, 0.012670, 6.265081e-06, 27.49, 1.2840, 23.55),
    (900.6, 1973.6576, 0.007500, 3.801475e-06, 16.27, 1.2931, 13.84),
    (950.5, 1949.1706, 0.006030, 3.092511e-06, 13.07, 1.2975, 11.08),
    (1000.0, 1924.8800, 0.004960, 2.575314e-06, 10.75, 1.3020, 9.08),
    (1200.0, 1826.7360, 0.002650, 1.451316e-06, 5.75, 1.3200, 4.79),
    # 1703 K: the reference mu value carries only 3 significant figures here, hence rtol 5e-3
    (1703.0, 1579.9038, 0.001050, 6.659128e-07, 2.28, 1.3653, 1.84),
]


@pytest.mark.parametrize("T,rho,mu,nu,pr_k11,k,pr_kT", TABLE)
def test_property_table(T, rho, mu, nu, pr_k11, k, pr_kT):
    f = FlibeFluid()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rtol_mu = 5e-3 if T > 1700.0 else 1e-3
        assert f.rho(T, 0) == pytest.approx(rho, rel=1e-4)
        assert f.mu(T, 0) == pytest.approx(mu, rel=rtol_mu)
        assert f.mu(T, 0) / f.rho(T, 0) == pytest.approx(nu, rel=rtol_mu)
        assert f.k(T, 0) == pytest.approx(k, rel=1e-4)
        assert f.Pr(T, 0) == pytest.approx(pr_kT, abs=0.02)
        # historical provenance column: Pr computed with a CONSTANT k = 1.1
        assert f.mu(T, 0) * f.cp(T, 0) / 1.1 == pytest.approx(pr_k11, abs=0.02)


def test_cp_is_constant():
    f = FlibeFluid()
    assert f.cp(800.0, 0) == 2386.0
    assert f.cp(1500.0, 0) == 2386.0


def test_hard_bounds_raise():
    f = FlibeFluid()
    for T in (730.0, 731.9, 1703.1, 2000.0):
        with pytest.raises(ValueError):
            f.rho(T, 0)
    f.rho(732.0, 0)          # inclusive at both ends
    f.rho(1703.0, 0)


def test_extrapolation_warns_inside_liquid_range():
    f = FlibeFluid()
    with pytest.warns(UserWarning, match="mu correlation extrapolated"):
        f.mu(1500.0, 0)
    with pytest.warns(UserWarning, match="k correlation extrapolated"):
        f.k(1500.0, 0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")       # inside the source range: must NOT warn
        f.mu(950.0, 0)
        f.k(1000.0, 0)


def test_viscosity_inversion():
    f = FlibeFluid()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert f.T_from_nu(3.09278e-6) == pytest.approx(950.478, abs=0.01)
