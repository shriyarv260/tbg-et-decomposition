import numpy as np
import pytest
from tbg_et.kinetics import log_rate, thermal_dos, KB_EV_K, trapz_weights


@pytest.fixture
def spectrum():
    e = np.linspace(-3, 3, 6001)
    return e, 0.1 + abs(e) + 4 * np.exp(-((e - 0.06)/0.03)**2)


@pytest.mark.parametrize("eta", [-0.2, 0, 0.15])
@pytest.mark.parametrize("mu", [-0.1, 0.12])
def test_pointwise_detailed_balance_for_asymmetric_dos(spectrum, eta, mu):
    e, d = spectrum
    red = log_rate(e, d, 0.8, 0.002, eta_v=eta, mu_ev=mu)
    ox = log_rate(e, d, 0.8, 0.002, eta_v=eta, mu_ev=mu, direction="oxidation")
    assert red - ox == pytest.approx(-eta / (KB_EV_K * 298.15), abs=1e-11)


def test_exact_coupling_and_dos_scaling(spectrum):
    e, d = spectrum
    base = log_rate(e, d, 0.8, 0.002)
    assert log_rate(e, d, 0.8, 0.006) - base == pytest.approx(np.log(9))
    assert log_rate(e, d * 7, 0.8, 0.002) - base == pytest.approx(np.log(7))
    assert log_rate(e, d, 0.8, 0.002, area_nm2=3) - base == pytest.approx(np.log(3))


def test_nonuniform_quadrature_weights():
    e = np.array([-2.0, -0.5, 0.2, 3.0])
    assert sum(trapz_weights(e)) == pytest.approx(5)
    assert np.sum(trapz_weights(e) * (2*e + 3)) == pytest.approx(20)


def test_energy_reference_shift_invariance(spectrum):
    e, d = spectrum
    assert log_rate(e, d, 1, 0.001, mu_ev=0.1) == pytest.approx(log_rate(e+2, d, 1, 0.001, mu_ev=2.1))


def test_vectorized_rates_match_scalar(spectrum):
    e, d = spectrum
    lam, h = np.array([0.7, 1.1]), np.array([0.001, 0.002])
    actual = log_rate(e, np.stack([d, 3*d]), lam, h)
    expected = [log_rate(e, d, lam[0], h[0]), log_rate(e, 3*d, lam[1], h[1])]
    np.testing.assert_allclose(actual, expected)


def test_zero_channels_and_extreme_barrier(spectrum):
    e, d = spectrum
    assert log_rate(e, d, 0.8, 0) == -np.inf
    assert log_rate(e, 0*d, 0.8, 0.002) == -np.inf
    assert np.isfinite(log_rate(e, d, 100, 0.001))


def test_constant_dos_high_driving_limit():
    e = np.linspace(-5, 5, 10001)
    # At strongly cathodic eta the normalized Gaussian lies below EF: integral -> DOS.
    from tbg_et.kinetics import HBAR_EV_S
    k = np.exp(log_rate(e, np.ones_like(e)*2, 0.8, 0.001, eta_v=-2))
    expected = 2*np.pi/HBAR_EV_S * 0.001**2 * 2
    assert k == pytest.approx(expected, rel=1e-7)


def test_thermal_dos_constant(spectrum):
    e, _ = spectrum
    assert thermal_dos(e, np.ones_like(e)*3.2) == pytest.approx(3.2)


@pytest.mark.parametrize("kwargs", [{"lambda_ev": -1}, {"coupling_ev": -1}, {"temperature_k": 0},
                                     {"area_nm2": 0}, {"eta_v": np.nan}, {"direction": "wrong"}])
def test_invalid_parameters(spectrum, kwargs):
    e, d = spectrum
    options = {"lambda_ev": 1, "coupling_ev": 0.001, **kwargs}
    with pytest.raises(ValueError):
        log_rate(e, d, **options)


def test_quadrature_convergence(spectrum):
    e, d = spectrum
    fine = log_rate(e, d, 0.82, 0.001)
    medium = log_rate(e[::2], d[::2], 0.82, 0.001)
    coarse = log_rate(e[::4], d[::4], 0.82, 0.001)
    # |E| has a cusp: trapezoidal convergence is quadratic, not spectral.
    assert abs(fine-medium) < 1e-5
    assert 3.8 < abs((medium-coarse)/(fine-medium)) < 4.2
