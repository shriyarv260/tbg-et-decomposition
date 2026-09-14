"""Standalone scientific PNG/PDF figures; no remote assets or web server."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .decomposition import FACTORS
from .kinetics import log_rate

COLORS = {"AA": "#d76b36", "AB": "#285c7a", "SP": "#739d85"}


def plot_analysis(data, summary, logs, baseline, ref, target, root):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "savefig.facecolor": "#fafaf8"})
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained", facecolor="#fafaf8")
    fig.suptitle("TBG / local electron-transfer decomposition", fontsize=22, x=0.03, ha="left")
    fig.supxlabel(f"{data.metadata['kind'].upper()} INPUTS  |  {summary['target_site']} / {summary['reference_site']}  |  "
                  "Nonadiabatic MHC-DOS · counterfactual attribution", fontsize=11)
    xy = data.xy_nm
    for ax, title, values, label, cmap in [
        (axes[0, 0], "Local stacking", np.array([list(COLORS).index(r) for r in data.registry]), "", None),
        (axes[0, 1], "Reorganization energy", data.lambda_ev, "λ (eV)", "viridis"),
        (axes[0, 2], "Full rate enhancement", (logs - logs[ref]) / np.log(10), "log₁₀(k / k reference)", "magma"),
    ]:
        if cmap is None:
            for label_, color in COLORS.items():
                mask = data.registry == label_
                ax.scatter(*xy[mask].T, c=color, s=14, label=label_, rasterized=True)
            ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.03, 1.04))
        else:
            m = ax.scatter(*xy.T, c=values, s=14, cmap=cmap, rasterized=True)
            fig.colorbar(m, ax=ax, shrink=0.8, label=label)
        ax.set(title=title, xlabel="x (nm)", ylabel="y (nm)", aspect="equal")
    ids = [ref, target]
    sp = np.flatnonzero(data.registry == "SP")
    if len(sp):
        ids.append(int(sp[len(sp)//2]))
    ids = list(dict.fromkeys(ids))
    ax = axes[1, 0]
    for i in ids:
        ax.plot(data.energy_ev, data.dos[i], label=data.site_ids[i], linewidth=1.7)
    mu = summary["rate_options"]["mu_ev"]
    left, right = max(data.energy_ev[0], mu - 0.3), min(data.energy_ev[-1], mu + 0.3)
    if left >= right:
        left, right = data.energy_ev[0], data.energy_ev[-1]
    ax.set(xlim=(left, right), title="Input local DOS", xlabel="Energy (eV)", ylabel="States / eV / nm²")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 1]
    eta = np.linspace(-0.35, 0.35, 71)
    opts = summary["rate_options"].copy()
    for i in ids[:2]:
        full, frozen = [], []
        for voltage in eta:
            opts["eta_v"] = voltage
            full.append(float(log_rate(data.energy_ev, data.dos[i], data.lambda_ev[i], data.coupling_ev[i], **opts)))
            frozen.append(float(log_rate(data.energy_ev, data.dos[i], data.lambda_ev[ref], data.coupling_ev[ref], **opts)))
        line, = ax.plot(eta, np.array(full)/np.log(10), label=str(data.site_ids[i]))
        ax.plot(eta, np.array(frozen)/np.log(10), ls="--", color=line.get_color(), alpha=0.7)
    ax.set(title="Potential sweep", xlabel="η (V); positive favors oxidation", ylabel="log₁₀ effective rate (s⁻¹)")
    ax.legend(frameon=False, fontsize=8)
    ax.text(0.03, 0.04, "Dashed: reference λ and coupling", transform=ax.transAxes, fontsize=8)
    ax = axes[1, 2]
    vals = [summary["decomposition"]["shapley_log_factors"][f] / np.log(10) for f in FACTORS]
    ax.bar(["DOS", "λ", "|H|²"], vals, color=["#285c7a", "#739d85", "#d76b36"], width=0.55)
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set(title="Attribution of enhancement", ylabel="Contribution to log₁₀ rate ratio")
    ax.text(0.97, 0.95, f"Total: {sum(vals):.3f} decades", transform=ax.transAxes, ha="right", va="top")
    for ax in axes.flat:
        ax.set_facecolor("#fafaf8")
    fig.savefig(root / "overview.png", dpi=160)
    fig.savefig(root / "overview.pdf")
    plt.close(fig)


def plot_twist(rows, output):
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    theta = [r["theta_deg"] for r in rows]
    for key, label, color in [("dos", "Direct DOS", "#285c7a"), ("lambda", "Reorganization", "#739d85"),
                              ("coupling", "Coupling", "#d76b36"), ("total", "Total", "#222222")]:
        ax.plot(theta, [r[key] for r in rows], "o-", label=label, color=color)
    ax.set(title="Synthetic AA/AB enhancement across twist angles", xlabel="Twist angle (degrees)",
           ylabel="Contribution to log₁₀(k AA / k AB)")
    ax.legend(frameon=False)
    fig.savefig(output, dpi=170)
    plt.close(fig)
