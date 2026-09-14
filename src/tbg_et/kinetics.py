"""Nonadiabatic, energy-independent-coupling MHC-DOS quadrature in eV.

Positive eta favors oxidation. Energy and mu use a common reference.
The returned rate is per redox molecule (s^-1), not a heterogeneous cm/s rate.
"""
import numpy as np

KB_EV_K = 8.617333262145e-5
HBAR_EV_S = 6.582119569e-16


def trapz_weights(energy):
    e = np.asarray(energy, dtype=float)
    if e.ndim != 1 or len(e) < 3 or not np.all(np.isfinite(e)) or np.any(np.diff(e) <= 0):
        raise ValueError("energy must contain at least 3 finite, strictly increasing values")
    widths = np.diff(e)
    return np.r_[widths[0] / 2, (widths[:-1] + widths[1:]) / 2, widths[-1] / 2]


def log_rate(energy, dos, lambda_ev, coupling_ev, *, temperature_k=298.15,
             eta_v=0.0, mu_ev=0.0, area_nm2=1.0, direction="reduction"):
    """Integrate areal DOS (states/eV/nm²) over an explicit coupling area.

    Leading dimensions of DOS broadcast with lambda/coupling; its last axis
    is energy. A zero coupling or identically zero DOS gives log(k)=-inf.
    Log-space integration avoids underflow of individual quadrature samples.
    """
    e = np.asarray(energy, dtype=float)
    weights = trapz_weights(e)
    d = np.asarray(dos, dtype=float)
    lam, h = np.broadcast_arrays(np.asarray(lambda_ev, float), np.asarray(coupling_ev, float))
    if d.ndim < 1 or d.shape[-1] != len(e) or not np.all(np.isfinite(d)) or np.any(d < 0):
        raise ValueError("DOS must be finite, nonnegative, with energy on the last axis")
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        raise ValueError("lambda_ev must be finite and positive")
    if not np.all(np.isfinite(h)) or np.any(h < 0):
        raise ValueError("coupling_ev must be a finite, nonnegative magnitude")
    if not np.isfinite(temperature_k) or temperature_k <= 0 or not np.isfinite(area_nm2) or area_nm2 <= 0:
        raise ValueError("temperature_k and area_nm2 must be finite and positive")
    if not np.isfinite(eta_v) or not np.isfinite(mu_ev):
        raise ValueError("eta_v and mu_ev must be finite")
    kt = KB_EV_K * temperature_k
    eps = e - mu_ev
    if direction == "reduction":
        occupancy = -np.logaddexp(0, eps / kt)
        gap = lam[..., None] + eta_v - eps
    elif direction == "oxidation":
        occupancy = -np.logaddexp(0, -eps / kt)
        gap = lam[..., None] - eta_v + eps
    else:
        raise ValueError("direction must be reduction or oxidation")
    with np.errstate(divide="ignore"):
        terms = (np.log(d) + np.log(weights) + occupancy
                 - gap**2 / (4 * lam[..., None] * kt))
        integral = np.logaddexp.reduce(terms, axis=-1)
        return (np.log(2 * np.pi / HBAR_EV_S) + 2 * np.log(h)
                + np.log(area_nm2) - 0.5 * np.log(4 * np.pi * lam * kt) + integral)


def thermal_dos(energy, dos, temperature_k=298.15, mu_ev=0.0):
    """Thermally broadened compressibility proxy, in states/eV/nm²."""
    weights = trapz_weights(energy)
    d = np.asarray(dos, float)
    if d.shape[-1] != len(energy) or np.any(d < 0) or not np.all(np.isfinite(d)):
        raise ValueError("invalid DOS")
    if not np.isfinite(temperature_k) or temperature_k <= 0 or not np.isfinite(mu_ev):
        raise ValueError("invalid temperature or chemical potential")
    kt = KB_EV_K * temperature_k
    x = (np.asarray(energy) - mu_ev) / kt
    derivative = np.exp(-np.logaddexp(0, x) - np.logaddexp(0, -x)) / kt
    return np.sum(d * derivative * weights, axis=-1)
