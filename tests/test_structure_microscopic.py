import numpy as np
import pytest
from tbg_et.structure import classify_registry, moire_period_nm
from tbg_et.microscopic import two_state_coupling, gap_reorganization, screening_surrogate


def test_reference_registries_and_periodicity():
    uv = np.array([[0, 0], [1/3, 1/3], [2/3, 2/3], [0.5, 0.5], [0, 0.5]])
    expected = ["AA", "AB", "AB", "SP", "SP"]
    assert list(classify_registry(uv)) == expected
    assert list(classify_registry(uv + [3, -4])) == expected
    assert moire_period_nm(1.1) == pytest.approx(12.8135, rel=1e-4)


def test_coupling_matches_matrix_lowdin_orthogonalization():
    h = np.array([[0.1, 0.035], [0.035, 0.4]])
    s = np.array([[1, 0.1], [0.1, 1]])
    values, vectors = np.linalg.eigh(s)
    invroot = (vectors * values**-0.5) @ vectors.T
    expected = abs((invroot @ h @ invroot)[0, 1])
    assert two_state_coupling(0.1, 0.4, 0.035, 0.1) == pytest.approx(expected)
    with pytest.raises(ValueError):
        two_state_coupling(0, 0, 0.1, 1)


def test_gap_estimate_and_bad_convention():
    result = gap_reorganization([0.9, 1.1], [-1.1, -0.9])
    assert result["lambda_mean_ev"] == pytest.approx(1)
    with pytest.raises(ValueError):
        gap_reorganization([-1, -0.9], [1, 0.9])


def test_surrogate_bounds_and_monotonicity():
    lam = screening_surrogate([0, 1, 10, 1e12])
    assert lam[0] == pytest.approx(1.17)
    assert lam[-1] == pytest.approx(0.82)
    assert np.all(np.diff(lam) < 0)
