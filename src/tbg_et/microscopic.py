"""Small extraction utilities for externally computed microscopic quantities."""
import numpy as np
from .kinetics import KB_EV_K


def two_state_coupling(h_aa_ev, h_bb_ev, h_ab_ev, overlap):
    """Löwdin-orthogonalized coupling for two real normalized diabatic states.

    H and S must belong to the SAME two-state basis. This is not a CP2K
    runner or a general projection-operator-diabatization implementation.
    """
    if not np.all(np.isfinite([h_aa_ev, h_bb_ev, h_ab_ev, overlap])) or abs(overlap) >= 1:
        raise ValueError("finite Hamiltonian elements and |overlap| < 1 required")
    return abs((h_ab_ev - overlap * (h_aa_ev + h_bb_ev) / 2) / (1 - overlap**2))


def gap_reorganization(gaps_oxidized_ev, gaps_reduced_ev, temperature_k=298.15):
    """Linear-response estimates with gap E_reduced - E_oxidized in BOTH ensembles.

    Report mean-gap and fluctuation estimates independently; disagreement is
    a diagnostic. Trajectory autocorrelation is not an uncertainty estimate.
    """
    arrays = [np.asarray(x, float) for x in (gaps_oxidized_ev, gaps_reduced_ev)]
    if any(x.ndim != 1 or len(x) < 2 or not np.all(np.isfinite(x)) for x in arrays):
        raise ValueError("each ensemble requires at least two finite energy gaps")
    if not np.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("positive finite temperature required")
    lam = (arrays[0].mean() - arrays[1].mean()) / 2
    if lam <= 0:
        raise ValueError("nonpositive mean-gap lambda: check ensemble labels and gap convention")
    return {"lambda_mean_ev": float(lam),
            "lambda_variance_oxidized_ev": float(arrays[0].var(ddof=1) / (2 * KB_EV_K * temperature_k)),
            "lambda_variance_reduced_ev": float(arrays[1].var(ddof=1) / (2 * KB_EV_K * temperature_k))}


def screening_surrogate(dos_thermal, *, lambda_metal_ev=0.82,
                        penalty_ev=0.35, dos_scale=1.0):
    """UNFITTED sensitivity model, not a reproduction of the 2024/2026 theory.

    lambda = lambda_metal + penalty / (1 + D_thermal / D_scale).
    All parameters must be declared in provenance; import calculated lambda
    for scientific inference. D_scale and D_thermal are states/eV/nm².
    """
    d = np.asarray(dos_thermal, float)
    if np.any(d < 0) or not np.all(np.isfinite(d)):
        raise ValueError("DOS proxy must be finite and nonnegative")
    if (not np.all(np.isfinite([lambda_metal_ev, penalty_ev, dos_scale]))
            or lambda_metal_ev <= 0 or penalty_ev < 0 or dos_scale <= 0):
        raise ValueError("invalid screening-surrogate parameters")
    return lambda_metal_ev + penalty_ev / (1 + d / dos_scale)
