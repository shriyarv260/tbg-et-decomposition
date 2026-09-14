import numpy as np
import pytest
from tbg_et.decomposition import decompose


def test_shapley_closure_and_coupling_identity():
    e = np.linspace(-2, 2, 2001)
    d = 0.1 + abs(e)
    target = d + 3 * np.exp(-(e/0.1)**2)
    result = decompose(e, (d, 1.2, 0.001), (target, 0.8, 0.002))
    assert result["closure_error"] == pytest.approx(0, abs=1e-12)
    assert result["shapley_log_factors"]["coupling"] == pytest.approx(np.log(4))
    assert len(result["hybrid_log_rates"]) == 8
    assert len(result["ordered_paths"]) == 6
    for path in result["ordered_paths"]:
        assert sum(path["log_factors"].values()) == pytest.approx(result["log_rate_ratio"])
    # The integral really is nonseparable for these inputs.
    assert np.ptp([p["log_factors"]["dos"] for p in result["ordered_paths"]]) > 0.01


def test_identical_inputs_and_reversal():
    e = np.linspace(-2, 2, 1601)
    a = (np.ones_like(e), 0.8, 0.001)
    b = (np.ones_like(e)*2, 1, 0.0015)
    assert decompose(e, a, a)["log_rate_ratio"] == 0
    forward, reverse = decompose(e, a, b), decompose(e, b, a)
    for factor in ("dos", "lambda", "coupling"):
        assert forward["shapley_log_factors"][factor] == pytest.approx(-reverse["shapley_log_factors"][factor])


def test_zero_reference_rejected():
    e = np.linspace(-2, 2, 1001)
    with pytest.raises(ValueError, match="positive finite"):
        decompose(e, (np.ones_like(e), 1, 0), (np.ones_like(e), 1, 0.001))
