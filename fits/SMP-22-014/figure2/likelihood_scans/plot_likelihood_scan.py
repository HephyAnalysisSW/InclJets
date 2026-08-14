#!/usr/bin/env python3
"""Plot direct and full-POD HERA/CMS likelihood scans from stored NPZ results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


GROUPS = (
    ("HERA", "HERA"),
    ("CMS", "CMS inclusive jets"),
    ("sum", "HERA + CMS"),
)


def make_figure(payload: dict[str, np.ndarray], parameter: str) -> plt.Figure:
    mask = payload["parameter_name"].astype(str) == parameter
    order = np.argsort(payload["scan_coordinate_sigma"][mask])
    x = payload["scan_coordinate_sigma"][mask][order]
    if x.size < 1:
        raise ValueError(f"No points stored for {parameter}")
    zero = int(np.argmin(np.abs(x)))

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(13.2, 7.0),
        sharex="col",
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        constrained_layout=True,
    )
    for column, (key, title) in enumerate(GROUPS):
        direct = payload[f"direct_chi2_{key}"][mask][order]
        pod = payload[f"full_pod_chi2_{key}"][mask][order]
        baseline = direct[zero]
        top = axes[0, column]
        bottom = axes[1, column]
        top.plot(x, direct - baseline, "o-", lw=1.8, ms=4.5, label="direct HERAPDF")
        top.plot(x, pod - baseline, "s--", lw=1.7, ms=4.2, label="full 100-mode POD")
        top.axhline(0.0, color="0.72", lw=0.8)
        top.axvline(0.0, color="0.82", lw=0.8)
        top.set_title(title)
        top.set_ylabel(r"$\chi^2-\chi^2_{\mathrm{direct}}(0)$")
        top.grid(alpha=0.2)
        delta = pod - direct
        bottom.plot(x, delta, "D-", color="#7a3e9d", lw=1.6, ms=4)
        bottom.axhline(0.0, color="0.35", lw=0.9)
        bottom.axvline(0.0, color="0.82", lw=0.8)
        bottom.set_xlabel(r"scan coordinate [$\sigma_{\mathrm{HESSE}}$]")
        bottom.set_ylabel(r"$\chi^2_{\rm POD}-\chi^2_{\rm direct}$")
        bottom.grid(alpha=0.2)
        bottom.text(
            0.04,
            0.92,
            rf"max $|\Delta\chi^2|={np.max(np.abs(delta)):.3g}$",
            transform=bottom.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
    axes[0, 0].legend(frameon=False, loc="best")
    values = payload["parameter_value"][mask][order]
    figure.suptitle(
        f"Unprofiled fixed-nuisance likelihood: {parameter}\n"
        f"{x.size} points; physical range {values.min():.7g} to {values.max():.7g}",
        fontsize=14,
    )
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with np.load(args.input, allow_pickle=False) as source:
        payload = {key: np.asarray(source[key]) for key in source.files}
    parameters = list(dict.fromkeys(payload["parameter_name"].astype(str).tolist()))
    pdf_path = output / "likelihood_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        for parameter in parameters:
            figure = make_figure(payload, parameter)
            figure.savefig(output / f"likelihood_{parameter}.png", dpi=180)
            pdf.savefig(figure)
            plt.close(figure)
    print(f"Wrote {len(parameters)} PNG(s) and {pdf_path}")


if __name__ == "__main__":
    main()
