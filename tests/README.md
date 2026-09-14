# Verification

Run `python -m pytest -q` from the repository root after installing `.[dev]`.

The tests cover detailed balance, quadrature convergence, exact coupling scaling,
decomposition closure, registry periodicity, microscopic extraction, input
validation, CSV joins, CLI execution, and figure generation. They verify the
software and its mathematical invariants, not the accuracy of synthetic physics.
