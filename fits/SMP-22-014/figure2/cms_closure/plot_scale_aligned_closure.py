#!/usr/bin/env python3
"""Make the scale-aligned POD PDF, jet-bin, and chi2 closure figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from plot_closure import pdf_figure


RAPIDITY = (r"$|y|<0.5$", r"$0.5<|y|<1.0$", r"$1.0<|y|<1.5$", r"$1.5<|y|<2.0$")


def cms_figure(path: Path, output: Path) -> None:
    with np.load(path, allow_pickle=False) as data:
        direct = {key.removeprefix("direct_input_"): data[key] for key in data.files if key.startswith("direct_input_")}
        pod = {key.removeprefix("all_flavours_"): data[key] for key in data.files if key.startswith("all_flavours_")}
    figure, axes = plt.subplots(2, 4, figsize=(15, 6.6), sharex="col", gridspec_kw={"height_ratios": [2.2, 1]}, constrained_layout=True)
    for index, (top, bottom) in enumerate(zip(axes[0], axes[1]), start=1):
        mask = direct["rapidity_bin"] == index
        pt, data_value, uncor = direct["pt"][mask], direct["data"][mask], direct["uncor"][mask]
        direct_theory, pod_theory = direct["theory"][mask], pod["theory"][mask]
        top.errorbar(pt, data_value, yerr=uncor, fmt="o", color="black", ms=3, lw=.8, label="CMS data (uncor.)")
        top.plot(pt, direct_theory, "o-", color="#1f77b4", ms=3, label="matched direct")
        top.plot(pt, pod_theory, "s--", color="#d62728", ms=3, label="full POD")
        top.set_xscale("log"); top.set_yscale("log"); top.set_title(RAPIDITY[index - 1]); top.grid(alpha=.2)
        bottom.plot(pt, 100 * (pod_theory / direct_theory - 1), "D-", color="#6a3d9a", ms=3)
        bottom.axhline(0, color=".35", lw=.8); bottom.set_xscale("log"); bottom.grid(alpha=.2)
        bottom.set_xlabel(r"$p_T$ [GeV]"); bottom.set_ylabel("POD/direct−1 [%]")
    axes[0, 0].set_ylabel(r"$d^2\sigma/dp_Tdy$"); axes[0, 0].legend(frameon=False, fontsize=8)
    figure.suptitle(r"CMS inclusive jets: scale-aligned direct versus full POD, $Q_0=1.65$ GeV", fontsize=14)
    figure.savefig(output / "scale_aligned_cms_predictions.png", dpi=180)
    plt.close(figure)


def chi2_figure(technical_path: Path, closure_path: Path, output: Path) -> None:
    technical, closure = yaml.safe_load(technical_path.read_text()), yaml.safe_load(closure_path.read_text())
    values = [
        technical["likelihoods"]["native"]["total_chi2"],
        technical["likelihoods"]["table_native_q0"]["total_chi2"],
        technical["likelihoods"]["table_q165"]["total_chi2"],
        closure["hybrid_likelihoods"]["all_flavours"]["total_chi2"],
    ]
    labels = ["native\ndirect", "direct table\n$Q_0=1.378$", "matched direct\n$Q_0=1.65$", "full POD\n$Q_0=1.65$"]
    increments = np.diff(values)
    figure, axis = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    baseline = values[0]
    axis.bar(np.arange(4), np.asarray(values) - baseline, bottom=baseline, color=["#777777", "#1f77b4", "#4c78a8", "#d62728"], width=.72)
    axis.axhline(baseline, color=".25", lw=.8)
    for index, value in enumerate(values):
        axis.text(index, value + .55, f"{value:.2f}", ha="center", va="bottom", fontsize=9)
    for index, increment in enumerate(increments, start=1):
        axis.annotate(f"{increment:+.2f}", xy=(index-.5, (values[index-1]+values[index])/2), ha="center", va="center", fontsize=9, bbox={"fc":"white", "ec":"none", "alpha":.85})
    axis.set_xticks(np.arange(4), labels); axis.set_ylabel(r"CMS-only fixed-zero-nuisance $\chi^2$")
    axis.set_title("Technical route and scale-aligned POD closure")
    axis.set_ylim(baseline - 2, values[-1] + 8)
    figure.savefig(output / "scale_aligned_chi2_waterfall.png", dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", required=True, type=Path, help="full_pod_projection.npz")
    parser.add_argument("--closure", required=True, type=Path, help="scale-aligned summary.yaml")
    parser.add_argument("--closure-data", required=True, type=Path, help="scale-aligned flavour_isolation.npz")
    parser.add_argument("--technical", required=True, type=Path, help="technical-controls summary.yaml")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdf_figure(args.projection, args.output_dir)
    cms_figure(args.closure_data, args.output_dir)
    chi2_figure(args.technical, args.closure, args.output_dir)
    print(f"Wrote scale-aligned closure figures in {args.output_dir}")


if __name__ == "__main__":
    main()
