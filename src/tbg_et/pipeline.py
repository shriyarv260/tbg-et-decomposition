"""End-to-end analysis with reproducible machine-readable artifacts."""
import csv
import hashlib
import json
from pathlib import Path
import platform
from importlib.metadata import version
import numpy as np
from . import __version__
from .kinetics import log_rate, thermal_dos
from .decomposition import decompose


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def site_index(data, site_id):
    hits = np.flatnonzero(data.site_ids == site_id)
    if len(hits) != 1:
        raise ValueError(f"unknown site ID: {site_id}")
    return int(hits[0])


def analyze(data, output, reference_id, target_id, *, temperature_k=298.15,
            eta_v=0.0, mu_ev=0.0, area_nm2=1.0, direction="reduction", plots=True):
    data.validate()
    ref, target = site_index(data, reference_id), site_index(data, target_id)
    options = dict(temperature_k=temperature_k, eta_v=eta_v, mu_ev=mu_ev,
                   area_nm2=area_nm2, direction=direction)
    e, d, lam, h = data.energy_ev, data.dos, data.lambda_ev, data.coupling_ev
    logs = log_rate(e, d, lam, h, **options)
    baseline = log_rate(e, d, lam[ref], h[ref], **options)
    decomposition = decompose(e, (d[ref], lam[ref], h[ref]), (d[target], lam[target], h[target]), **options)
    # Coarsening tests quadrature sensitivity only, not energy-window sufficiency.
    indices = np.unique(np.r_[np.arange(0, len(e), 2), len(e) - 1])
    coarse = log_rate(e[indices], d[:, indices], lam, h, **options)
    if not all(np.all(np.isfinite(x)) for x in (logs, baseline, coarse)):
        raise ValueError("nonfinite rates: check energy window, positive channels, and parameter scales")
    coarse_delta = float(np.max(np.abs(coarse - logs)))
    warnings = ["Numerical checks do not establish accuracy of physical input models.",
                "Verify energy-window coverage and upstream spectral resolution independently.",
                "Rate units are effective per-molecule s^-1; no conversion to SECCM current or cm/s."]
    if data.metadata["kind"] == "synthetic":
        warnings.insert(0, "SYNTHETIC DEMONSTRATION: no experimental or first-principles predictions.")
    if coarse_delta > 0.01:
        warnings.append("Coarsened-grid log-rate change exceeds 0.01; refine the energy grid.")
    if not e[0] < mu_ev < e[-1]:
        warnings.append("Chemical potential lies outside the supplied energy window.")
    root = Path(output)
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"output directory is not empty: {root}; choose a new directory")
    root.mkdir(parents=True, exist_ok=True)
    data.save(root)
    summary = {
        "version": __version__, "input_kind": data.metadata["kind"],
        "reference_site": reference_id, "target_site": target_id,
        "n_sites": len(data.site_ids), "rate_options": options,
        "reference_log_rate_s-1": float(logs[ref]), "target_log_rate_s-1": float(logs[target]),
        "dos_only_log_ratio": float(baseline[target] - baseline[ref]),
        "decomposition": decomposition,
        "numerics": {"coarsening_max_abs_log_rate_change": coarse_delta,
                     "energy_min_ev": float(e[0]), "energy_max_ev": float(e[-1]),
                     "max_energy_step_ev": float(np.max(np.diff(e)))},
        "warnings": warnings,
    }
    write_json(root / "summary.json", summary)
    thermo = thermal_dos(e, d, temperature_k, mu_ev)
    with (root / "rates.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["site_id", "registry", "x_nm", "y_nm", "lambda_ev", "coupling_ev",
                         "thermal_dos_states_per_ev_nm2", "log_rate_s-1", "log10_rate_ratio_to_reference",
                         "dos_only_log10_ratio_to_reference"])
        for i, name in enumerate(data.site_ids):
            writer.writerow([name, data.registry[i], *data.xy_nm[i], lam[i], h[i], thermo[i], logs[i],
                             (logs[i] - logs[ref]) / np.log(10), (baseline[i] - baseline[ref]) / np.log(10)])
    if plots:
        from .plotting import plot_analysis
        plot_analysis(data, summary, logs, baseline, ref, target, root)
    attribution = decomposition["shapley_log_factors"]
    text = f"""# Local electron-transfer decomposition

Input status: **{data.metadata['kind'].upper()}**. {warnings[0]}

Comparison: `{target_id}` / `{reference_id}` at {temperature_k:g} K,
eta = {eta_v:g} V, mu = {mu_ev:g} eV, effective coupling area = {area_nm2:g} nm².

| Contribution | log10 enhancement |
| --- | ---: |
| Direct DOS | {attribution['dos'] / np.log(10):.6f} |
| Reorganization energy | {attribution['lambda'] / np.log(10):.6f} |
| Coupling | {attribution['coupling'] / np.log(10):.6f} |
| Total | {decomposition['log10_rate_ratio']:.6f} |

Values use the average over all six counterfactual orders in log-rate space.
They are descriptive attributions, not independent causal mechanisms. DOS is
held separate from lambda during each intervention. All paths and eight hybrid
rates are in `summary.json`; maps and spectra are in `overview.png` when plotted.

Closure error: {decomposition['closure_error']:.3e} (natural log units).
Maximum log-rate change on a coarsened energy grid: {coarse_delta:.3e}.

## Limitations

""" + "\n".join(f"- {w}" for w in warnings + data.metadata.get("limitations", [])) + "\n"
    (root / "report.md").write_text(text)
    update_manifest(root)
    return summary


def update_manifest(root):
    root = Path(root)
    files = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(root.rglob("*")) if p.is_file() and p != root / "manifest.json"}
    write_json(root / "manifest.json", {"package_version": __version__,
               "python": platform.python_version(), "numpy": np.__version__,
               "matplotlib": version("matplotlib"), "sha256": files})
