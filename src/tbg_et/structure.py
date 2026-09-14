"""Registry classification from an externally supplied local displacement field."""
import numpy as np

# Fractional coordinates in graphene's 60-degree primitive basis.
REFERENCES = {"AA": [(0, 0)], "AB": [(1/3, 1/3), (2/3, 2/3)],
              "SP": [(0.5, 0), (0, 0.5), (0.5, 0.5)]}


def moire_period_nm(theta_deg, lattice_nm=0.246):
    if not np.isfinite(theta_deg) or not 0 < theta_deg <= 30:
        raise ValueError("twist angle must lie in (0, 30] degrees")
    return lattice_nm / (2 * np.sin(np.deg2rad(theta_deg) / 2))


def registry_distances(displacement_fractional):
    uv = np.asarray(displacement_fractional, float)
    if uv.ndim != 2 or uv.shape[1] != 2 or not np.all(np.isfinite(uv)):
        raise ValueError("displacement must be a finite (n_sites, 2) array")
    basis = np.array([[1, 0], [0.5, np.sqrt(3) / 2]])
    uv = uv % 1
    distances = []
    for refs in REFERENCES.values():
        images = np.array([(a + i, b + j) for a, b in refs
                           for i in (-1, 0, 1) for j in (-1, 0, 1)])
        delta = (uv[:, None, :] - images[None, :, :]) @ basis
        distances.append(np.sqrt(np.min(np.sum(delta**2, axis=-1), axis=1)))
    return np.stack(distances, axis=-1)


def classify_registry(displacement_fractional):
    """Nearest ideal registry; AB includes BA. Not an atomistic relaxation."""
    return np.array(list(REFERENCES))[registry_distances(displacement_fractional).argmin(axis=1)]
