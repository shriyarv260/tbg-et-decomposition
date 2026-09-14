"""Counterfactual factor replacement; descriptive attribution, not causality."""
from itertools import permutations
import numpy as np
from .kinetics import log_rate

FACTORS = ("dos", "lambda", "coupling")


def decompose(energy, reference, target, **rate_options):
    """Evaluate all eight AB/target hybrids and exact log-rate Shapley values.

    reference and target are (dos, lambda_ev, coupling_ev). Lambda is frozen
    during a DOS intervention even when the supplied lambda came from DOS.
    This separates direct DOS and mediated lambda effects by convention.
    """
    logs = {}
    for mask in range(8):
        args = [target[i] if mask & (1 << i) else reference[i] for i in range(3)]
        value = float(log_rate(energy, *args, **rate_options))
        if not np.isfinite(value):
            raise ValueError("decomposition requires strictly positive finite rates in all eight hybrids")
        logs[mask] = value
    paths = []
    contributions = np.zeros(3)
    for order in permutations(range(3)):
        mask = 0
        steps = {}
        for i in order:
            next_mask = mask | (1 << i)
            delta = logs[next_mask] - logs[mask]
            contributions[i] += delta / 6
            steps[FACTORS[i]] = delta
            mask = next_mask
        paths.append({"order": [FACTORS[i] for i in order], "log_factors": steps})
    total = logs[7] - logs[0]
    return {
        "log_rate_ratio": total,
        "log10_rate_ratio": total / np.log(10),
        "shapley_log_factors": dict(zip(FACTORS, contributions.tolist())),
        "closure_error": float(sum(contributions) - total),
        "hybrid_log_rates": {str(k): v for k, v in logs.items()},
        "ordered_paths": paths,
    }
