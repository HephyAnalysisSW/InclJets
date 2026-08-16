#!/usr/bin/env python3
"""Make PDF and CMS-jet closure figures from stored POD diagnostic arrays."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FLAVOR_LABELS = {21: "g", 2: "u", -2: r"$\bar u$", 1: "d", -1: r"$\bar d$", 3: "s", -3: r"$\bar s$", 4: "c", -4: r"$\bar c$", 5: "b", -5: r"$\bar b$"}
RAPIDITY_LABELS = (r"$|y|<0.5$", r"$0.5<|y|<1.0$", r"$1.0<|y|<1.5$", r"$1.5<|y|<2.0$")


def pdf_figure(path: Path, output: Path) -> None:
    with np.load(path, allow_pickle=False) as data:
        x, flavors = data["x_grid"], data["flavors"].astype(int)
        target, projected = data["target_grid"], data["projected_grid"]
    figure, axes = plt.subplots(3, 4, figsize=(13, 8), sharex=True)
    global_scale = max(float(np.max(np.abs(target))), 1e-30)
    for axis, pid, direct, pod in zip(axes.flat, flavors, target, projected):
        scale = float(np.max(np.abs(direct)))
        if scale < 1e-12 * global_scale:
            axis.set_title(FLAVOR_LABELS.get(int(pid), str(pid)))
            axis.text(0.5, 0.5, "inactive at this scale", ha="center", va="center", transform=axis.transAxes, fontsize=9)
            axis.set_xscale("log")
            continue
        # Ratios are ill-defined wherever an individual flavour vanishes or
        # changes sign at x -> 1.  A per-flavour normalised difference keeps
        # every curve finite and directly ranks the relevant closure defect.
        difference = (pod - direct) / scale
        axis.semilogx(x, difference, color="#6a3d9a", lw=1.4)
        axis.axhline(0, color="0.35", lw=0.8)
        axis.set_title(FLAVOR_LABELS.get(int(pid), str(pid)))
        axis.grid(alpha=0.2)
        max_percent = 100 * np.max(np.abs(difference))
        axis.text(0.04, 0.08, f"max {max_percent:.2g}%", transform=axis.transAxes, fontsize=8)
    axes.flat[-1].set_visible(False)
    for axis in axes[-1, :]: axis.set_xlabel(r"$x$")
    for axis in axes[:, 0]: axis.set_ylabel(r"$(f_{\rm POD}-f_{\rm direct})/\max_x|f_{\rm direct}|$")
    figure.suptitle(r"Full-100-mode POD PDF closure at $Q=1.65$ GeV", fontsize=14)
    figure.tight_layout()
    figure.savefig(output / "pdf_closure_by_flavor.png", dpi=180)
    plt.close(figure)


def cms_figures(path: Path, output: Path) -> None:
    with np.load(path, allow_pickle=False) as data:
        direct = {key.removeprefix("direct_"): np.asarray(data[key]) for key in data.files if key.startswith("direct_") and key not in {"direct_total_chi2", "direct_data_chi2", "direct_log_penalty_chi2"}}
        pod = {key.removeprefix("pod_"): np.asarray(data[key]) for key in data.files if key.startswith("pod_") and key not in {"pod_total_chi2", "pod_data_chi2", "pod_log_penalty_chi2"}}
        chi2 = {key: float(data[key]) for key in data.files if key.endswith("chi2")}
    figure, axes = plt.subplots(2, 4, figsize=(15, 6.6), sharex="col", gridspec_kw={"height_ratios": [2.2, 1]}, constrained_layout=True)
    residual_figure, residual_axes = plt.subplots(1, 4, figsize=(15, 3.6), sharex=True, constrained_layout=True)
    for index, (top, bottom, residual) in enumerate(zip(axes[0], axes[1], residual_axes)):
        mask = direct["rapidity_bin"] == index + 1
        pt, data_value, uncertainty = direct["pt"][mask], direct["data"][mask], direct["uncor"][mask]
        direct_theory, pod_theory = direct["theory"][mask], pod["theory"][mask]
        top.errorbar(pt, data_value, yerr=uncertainty, fmt="o", color="black", ms=3, lw=0.8, label="CMS data (uncor.)")
        top.plot(pt, direct_theory, "o-", color="#1f77b4", ms=3, label="direct HERAPDF")
        top.plot(pt, pod_theory, "s--", color="#d62728", ms=3, label="full POD")
        top.set_xscale("log"); top.set_yscale("log"); top.set_title(RAPIDITY_LABELS[index]); top.grid(alpha=0.2)
        bottom.plot(pt, pod_theory / direct_theory - 1, "D-", color="#6a3d9a", ms=3)
        bottom.axhline(0, color="0.35", lw=0.8); bottom.set_xscale("log"); bottom.grid(alpha=0.2)
        bottom.set_xlabel(r"$p_T$ [GeV]"); bottom.set_ylabel(r"POD/direct$-1$")
        delta_pull2 = ((data_value - pod_theory) / uncertainty) ** 2 - ((data_value - direct_theory) / uncertainty) ** 2
        residual.axhline(0, color="0.35", lw=0.8); residual.plot(pt, delta_pull2, "D-", color="#6a3d9a", ms=3)
        residual.set_xscale("log"); residual.grid(alpha=0.2); residual.set_title(RAPIDITY_LABELS[index]); residual.set_xlabel(r"$p_T$ [GeV]")
    axes[0, 0].set_ylabel(r"$d^2\sigma/dp_Tdy$"); axes[0, 0].legend(frameon=False, fontsize=8)
    figure.suptitle("CMS inclusive jets: direct versus POD, all nuisance shifts fixed to zero", fontsize=14)
    figure.savefig(output / "cms_zero_systematics_predictions.png", dpi=180)
    plt.close(figure)
    residual_axes[0].set_ylabel(r"$\Delta[(D-T)^2/\sigma_{\rm uncor}^2]$")
    residual_figure.suptitle("CMS bin-level change in diagonal residual proxy (systematics fixed to zero)", fontsize=13)
    residual_figure.savefig(output / "cms_zero_systematics_residual_proxy.png", dpi=180)
    plt.close(residual_figure)
    (output / "figure_summary.txt").write_text(
        "CMS-only, fixed-zero-nuisance likelihoods\n"
        f"direct total/data/log = {chi2['direct_total_chi2']:.12g} / {chi2['direct_data_chi2']:.12g} / {chi2['direct_log_penalty_chi2']:.12g}\n"
        f"POD    total/data/log = {chi2['pod_total_chi2']:.12g} / {chi2['pod_data_chi2']:.12g} / {chi2['pod_log_penalty_chi2']:.12g}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", required=True, type=Path, help="full_pod_projection.npz")
    parser.add_argument("--cms", required=True, type=Path, help="cms_zero_systematics.npz")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdf_figure(args.projection, args.output_dir)
    cms_figures(args.cms, args.output_dir)
    print(f"Wrote closure figures in {args.output_dir}")


if __name__ == "__main__":
    main()
