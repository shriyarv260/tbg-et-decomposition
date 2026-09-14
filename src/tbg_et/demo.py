"""Deterministic synthetic fixture. No electronic-structure solver is run."""
import numpy as np
from .data import Dataset
from .kinetics import thermal_dos
from .microscopic import screening_surrogate
from .structure import classify_registry, registry_distances, moire_period_nm


def make_demo(theta_deg=1.1, grid_size=24, energy_points=1601, temperature_k=298.15):
    if not isinstance(grid_size, int) or grid_size < 6 or grid_size % 6:
        raise ValueError("grid_size must be a positive multiple of 6, at least 6")
    if not isinstance(energy_points, int) or energy_points < 401 or energy_points % 2 != 1:
        raise ValueError("energy_points must be odd and at least 401")
    period = moire_period_nm(theta_deg)
    u, v = np.meshgrid(np.arange(grid_size) / grid_size, np.arange(grid_size) / grid_size)
    uv = np.column_stack((u.ravel(), v.ravel()))
    xy = period * uv @ np.array([[1, 0], [0.5, np.sqrt(3)/2]])
    labels = classify_registry(uv)
    aa = np.exp(-(registry_distances(uv)[:, 0] / 0.14)**2)
    magic = np.exp(-((theta_deg - 1.1) / 0.45)**2)
    energy = np.linspace(-2, 2, energy_points)
    # Synthetic positive areal DOS, not a continuum or tight-binding spectrum.
    background = 0.15 + 1.4 * np.abs(energy)
    peak = np.exp(-0.5 * ((energy - 0.045) / 0.028)**2) + np.exp(-0.5 * ((energy + 0.045) / 0.028)**2)
    dos = background[None, :] + (0.12 + 3.0 * aa[:, None] * magic) * peak[None, :]
    dthermal = thermal_dos(energy, dos, temperature_k)
    lam = screening_surrogate(dthermal)
    coupling = 0.001 * (1 + 0.35 * aa * magic)
    ids = np.array([f"site-{i:04d}" for i in range(len(uv))])
    ids[0] = "AA-center"
    ids[(grid_size // 3) * grid_size + grid_size // 3] = "AB-center"
    ids[(grid_size // 2) * grid_size + grid_size // 2] = "SP-center"
    return Dataset(energy, dos, xy, labels, lam, coupling, ids, {
        "schema_version": 1, "kind": "synthetic", "source": "tbg_et.demo.make_demo v0.1.0",
        "dos_units": "states/eV/nm^2", "energy_reference": "synthetic charge-neutrality point at 0 eV",
        "structure_model": "rigid periodic displacement grid; NOT relaxed atomic coordinates",
        "dos_model": "positive linear background plus two Gaussian peaks; NOT computed LDOS",
        "lambda_model": "unfitted bounded screening surrogate; NOT the published nonlocal dielectric model",
        "coupling_model": "prescribed 1 meV baseline with synthetic AA modulation; NOT CDFT/POD",
        "theta_deg": theta_deg, "grid_size": grid_size, "energy_points": energy_points,
        "temperature_k_for_screening": temperature_k, "moire_period_nm": period,
        "parameters": {"lambda_metal_ev": 0.82, "screening_penalty_ev": 0.35,
                       "dos_scale_states_per_ev_nm2": 1.0, "aa_width_fractional": 0.14,
                       "magic_center_deg": 1.1, "magic_width_deg": 0.45},
        "limitations": ["No relaxed geometry, computed wavefunctions, or measured rate data",
                        "Synthetic twist enhancement is prescribed, not predicted",
                        "No validation of the AA anomaly or claim of research novelty"],
    }).validate()
