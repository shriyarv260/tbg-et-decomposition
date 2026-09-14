"""Validated, pickle-free interchange format for independent upstream solvers."""
from dataclasses import dataclass
import csv
import json
from pathlib import Path
import numpy as np
from .kinetics import trapz_weights

ARRAY_KEYS = ("energy_ev", "dos", "xy_nm", "registry", "lambda_ev", "coupling_ev", "site_ids")


@dataclass
class Dataset:
    energy_ev: np.ndarray
    dos: np.ndarray
    xy_nm: np.ndarray
    registry: np.ndarray
    lambda_ev: np.ndarray
    coupling_ev: np.ndarray
    site_ids: np.ndarray
    metadata: dict

    def validate(self):
        trapz_weights(self.energy_ev)
        n = len(self.site_ids)
        if n < 2 or self.site_ids.shape != (n,) or self.site_ids.dtype.kind not in "US":
            raise ValueError("site_ids must be a one-dimensional string array with at least two sites")
        if len(set(self.site_ids.tolist())) != n or any(not s.strip() for s in self.site_ids):
            raise ValueError("site_ids must be unique and nonempty")
        if self.dos.shape != (n, len(self.energy_ev)) or not np.all(np.isfinite(self.dos)) or np.any(self.dos < 0):
            raise ValueError("dos must be finite nonnegative (n_sites, n_energy)")
        if np.any(np.sum(self.dos, axis=1) == 0):
            raise ValueError("each site needs nonzero DOS for decomposition")
        if self.xy_nm.shape != (n, 2) or not np.all(np.isfinite(self.xy_nm)):
            raise ValueError("xy_nm must be finite (n_sites, 2)")
        if self.registry.shape != (n,) or not set(self.registry.tolist()) <= {"AA", "AB", "SP"}:
            raise ValueError("registry must contain only AA, AB (including BA), or SP")
        for key in ("lambda_ev", "coupling_ev"):
            x = getattr(self, key)
            if x.shape != (n,) or not np.all(np.isfinite(x)) or np.any(x <= 0):
                raise ValueError(f"{key} must be a positive finite (n_sites,) array")
        m = self.metadata
        if not isinstance(m, dict) or m.get("schema_version") != 1:
            raise ValueError("metadata schema_version must be 1")
        if m.get("kind") not in ("synthetic", "external"):
            raise ValueError("metadata kind must be synthetic or external")
        if m.get("dos_units") != "states/eV/nm^2":
            raise ValueError("DOS units must be states/eV/nm^2; convert upstream normalization explicitly")
        for key in ("source", "energy_reference", "structure_model", "dos_model", "lambda_model", "coupling_model"):
            if not isinstance(m.get(key), str) or not m[key].strip():
                raise ValueError(f"metadata requires a nonempty {key}")
        # Fail early on non-JSON provenance, including NaN/Infinity.
        json.dumps(m, allow_nan=False)
        return self

    def save(self, directory):
        self.validate()
        path = Path(directory)
        path.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path / "input.npz", **{k: getattr(self, k) for k in ARRAY_KEYS})
        (path / "metadata.json").write_text(json.dumps(self.metadata, indent=2, allow_nan=False) + "\n")


def load_dataset(npz_path, metadata_path):
    with np.load(npz_path, allow_pickle=False) as archive:
        if not set(ARRAY_KEYS) <= set(archive.files):
            raise ValueError(f"NPZ requires arrays: {', '.join(ARRAY_KEYS)}")
        arrays = {k: archive[k] for k in ARRAY_KEYS}
    return Dataset(**arrays, metadata=json.loads(Path(metadata_path).read_text())).validate()


def load_csv(sites_path, dos_path, metadata_path):
    """Join sites to long-form LDOS by IDs; reject missing/duplicate energies."""
    with Path(sites_path).open(newline="") as file:
        sites = list(csv.DictReader(file))
    if not sites:
        raise ValueError("sites CSV is empty")
    ids = [r["site_id"] for r in sites]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate site IDs in sites CSV")
    spectra = {name: {} for name in ids}
    with Path(dos_path).open(newline="") as file:
        for row in csv.DictReader(file):
            name, e, value = row["site_id"], float(row["energy_ev"]), float(row["dos_states_per_ev_nm2"])
            if name not in spectra:
                raise ValueError(f"LDOS references unknown site: {name}")
            if e in spectra[name]:
                raise ValueError(f"duplicate energy for site {name}: {e}")
            if not np.all(np.isfinite([e, value])):
                raise ValueError("nonfinite LDOS sample")
            spectra[name][e] = value
    energy = sorted(spectra[ids[0]])
    if any(sorted(spectra[name]) != energy for name in ids):
        raise ValueError("every site must have the same complete energy grid")
    return Dataset(np.array(energy), np.array([[spectra[name][e] for e in energy] for name in ids]),
                   np.array([[float(r["x_nm"]), float(r["y_nm"])] for r in sites]),
                   np.array([r["registry"] for r in sites]),
                   np.array([float(r["lambda_ev"]) for r in sites]),
                   np.array([float(r["coupling_ev"]) for r in sites]), np.array(ids),
                   json.loads(Path(metadata_path).read_text())).validate()
