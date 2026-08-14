#!/usr/bin/env python3
"""Plot an external Hessian PDF set and its saved POD projection."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


FLAVOR_LABELS = {
    -5: r"$\bar b$",
    -4: r"$\bar c$",
    -3: r"$\bar s$",
    -2: r"$\bar u$",
    -1: r"$\bar d$",
    1: r"$d$",
    2: r"$u$",
    3: r"$s$",
    4: r"$c$",
    5: r"$b$",
    21: r"$g$",
}


def parse_args() -> argparse.Namespace:
    here = Path(__file__).parent
    parser = argparse.ArgumentParser(
        description="Plot Hessian PDF bands before and after POD projection."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=here / "outputs" / "ct18nnlo_pod_projection.npz",
    )
    parser.add_argument("--output-directory", type=Path, default=here / "plots")
    return parser.parse_args()


def finite_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    out = np.full_like(numerator, np.nan, dtype=float)
    # Match the upstream demo: suppress ratios in tails where the denominator
    # is below 1e-8 of its flavor-wise maximum.
    scale = max(float(np.max(np.abs(denominator))), 1e-14)
    active = np.abs(denominator) > 1e-8 * scale
    out[active] = numerator[active] / denominator[active]
    return out


def main() -> None:
    args = parse_args()
    with np.load(args.input) as result:
        metadata = json.loads(str(result["metadata"]))
        x_grid = result["x_grid"]
        flavors = result["flavors"]
        target = result["target_grids"][0]
        projected = result["projected_grids"][0]
        target_error = result["target_hessian_band"]
        projected_error = result["projected_hessian_band"]
        relative_residuals = result["relative_residuals"]
        coefficient_covariance = result["coefficient_covariance"]

    if not target_error.size:
        raise ValueError("Saved projection does not contain a Hessian band")

    args.output_directory.mkdir(parents=True, exist_ok=True)
    n_columns = 3
    n_rows = int(np.ceil(len(flavors) / n_columns))

    fig, axes = plt.subplots(n_rows, n_columns, figsize=(13, 3.2 * n_rows))
    for index, (axis, pid) in enumerate(zip(axes.flat, flavors)):
        axis.fill_between(
            x_grid,
            target[index] - target_error[index],
            target[index] + target_error[index],
            color="tab:blue",
            alpha=0.25,
            label=f"{metadata['target_set']} Hessian",
        )
        axis.fill_between(
            x_grid,
            projected[index] - projected_error[index],
            projected[index] + projected_error[index],
            color="tab:orange",
            alpha=0.25,
            label="POD projection",
        )
        axis.plot(x_grid, target[index], color="tab:blue", lw=1.1)
        axis.plot(x_grid, projected[index], color="tab:orange", lw=1.1)
        axis.set_xscale("log")
        axis.set_xlim(x_grid[0], x_grid[-1])
        axis.set_title(FLAVOR_LABELS.get(int(pid), f"pid {pid}"))
        axis.set_xlabel("x")
        axis.set_ylabel(r"$x f(x,Q)$")
        axis.grid(alpha=0.2)
    for axis in axes.flat[len(flavors) :]:
        axis.set_visible(False)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(
        f"{metadata['target_set']} projected on {metadata['basis_set']} "
        f"at Q={metadata['q_GeV']:g} GeV"
    )
    fig.tight_layout()
    band_png = args.output_directory / "ct18nnlo_pod_bands.png"
    band_pdf = args.output_directory / "ct18nnlo_pod_bands.pdf"
    fig.savefig(band_png, dpi=180)
    fig.savefig(band_pdf)
    plt.close(fig)

    fig, axes = plt.subplots(n_rows, n_columns, figsize=(13, 3.2 * n_rows))
    for index, (axis, pid) in enumerate(zip(axes.flat, flavors)):
        central_ratio = finite_ratio(projected[index], target[index])
        error_ratio = finite_ratio(projected_error[index], target_error[index])
        axis.plot(x_grid, central_ratio, color="tab:green", label="central ratio")
        axis.plot(x_grid, error_ratio, color="tab:red", label="uncertainty ratio")
        axis.axhline(1.0, color="black", lw=0.7, ls="--")
        axis.set_xscale("log")
        # Some heavy-flavor grids vanish identically at Q=1.65 GeV. In that
        # case both ratios are undefined and Matplotlib cannot infer log limits.
        axis.set_xlim(x_grid[0], x_grid[-1])
        axis.set_title(FLAVOR_LABELS.get(int(pid), f"pid {pid}"))
        axis.set_xlabel("x")
        axis.set_ylabel("projected / target")
        finite = np.concatenate(
            (central_ratio[np.isfinite(central_ratio)], error_ratio[np.isfinite(error_ratio)])
        )
        if finite.size:
            y_min = max(0.0, float(np.min(finite)) - 0.05)
            y_max = min(2.0, float(np.max(finite)) + 0.05)
            if y_min >= y_max:
                y_min, y_max = 0.95, 1.05
        else:
            y_min, y_max = 0.95, 1.05
            axis.text(0.5, 0.5, "inactive at this Q", ha="center", transform=axis.transAxes)
        axis.set_ylim(y_min, y_max)
        axis.grid(alpha=0.2)
    for axis in axes.flat[len(flavors) :]:
        axis.set_visible(False)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle("POD reconstruction ratios (vanishing-PDF tails masked)")
    fig.tight_layout()
    ratio_png = args.output_directory / "ct18nnlo_pod_ratios.png"
    ratio_pdf = args.output_directory / "ct18nnlo_pod_ratios.pdf"
    fig.savefig(ratio_png, dpi=180)
    fig.savefig(ratio_pdf)
    plt.close(fig)

    # This is the upstream plot named ``ratio_bands``: both Hessian bands and
    # both central curves are normalized to the *CT18 central* member. It is
    # distinct from the central/error reconstruction diagnostic above.
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(13, 3.2 * n_rows))
    for index, (axis, pid) in enumerate(zip(axes.flat, flavors)):
        denominator = target[index]
        denominator_scale = max(float(np.max(np.abs(denominator))), 1e-14)
        active = np.abs(denominator) > 1e-8 * denominator_scale
        if np.any(active):
            x_active = x_grid[active]
            target_ratio = np.ones_like(x_active)
            target_ratio_error = target_error[index][active] / np.abs(
                denominator[active]
            )
            projected_ratio = projected[index][active] / denominator[active]
            projected_ratio_error = projected_error[index][active] / np.abs(
                denominator[active]
            )

            finite = np.concatenate(
                (
                    target_ratio - target_ratio_error,
                    target_ratio + target_ratio_error,
                    projected_ratio - projected_ratio_error,
                    projected_ratio + projected_ratio_error,
                )
            )
            finite = finite[np.isfinite(finite)]
            y_min = max(-1.0, float(np.min(finite)) - 0.05) if finite.size else 0.8
            y_max = min(3.0, float(np.max(finite)) + 0.05) if finite.size else 1.2
            if y_min >= y_max:
                y_min, y_max = 0.8, 1.2

            axis.fill_between(
                x_active,
                target_ratio - target_ratio_error,
                target_ratio + target_ratio_error,
                facecolor="0.75",
                edgecolor="0.55",
                linewidth=0.4,
            )
            axis.fill_between(
                x_active,
                projected_ratio - projected_ratio_error,
                projected_ratio + projected_ratio_error,
                facecolor="none",
                edgecolor="tab:orange",
                hatch="//////",
                linewidth=0.0,
            )
            axis.plot(x_active, target_ratio, color="black", lw=1.0)
            axis.plot(
                x_active, projected_ratio, color="tab:red", lw=1.0, ls="--"
            )
        else:
            y_min, y_max = 0.95, 1.05
            axis.text(
                0.5, 0.5, "inactive", ha="center", va="center", transform=axis.transAxes
            )

        axis.axhline(1.0, color="black", lw=0.6)
        axis.set_xscale("log")
        axis.set_xlim(x_grid[0], x_grid[-1])
        axis.set_ylim(y_min, y_max)
        axis.set_title(FLAVOR_LABELS.get(int(pid), f"pid {pid}"))
        axis.set_xlabel("x")
        axis.set_ylabel(f"PDF / {metadata['target_set']} central")
        axis.grid(alpha=0.15)

    covariance_rank = int(np.linalg.matrix_rank(coefficient_covariance))
    coefficient_sigma = np.sqrt(
        np.clip(np.diag(coefficient_covariance), 0.0, None)
    )
    info_axis = axes.flat[-1]
    info_axis.set_axis_off()
    info_axis.text(
        0.05,
        0.92,
        "\n".join(
            (
                f"{metadata['target_set']}: {metadata['hessian_pairs']} pairs, "
                f"{metadata['hessian_confidence_level']}% CL",
                f"{metadata['n_basis']} POD modes, {metadata['metric']}",
                f"central res/shift {relative_residuals[0]:.2e}",
                f"rank(Cc) {covariance_rank}/{metadata['n_basis']}",
                f"rms/max sigma(c) "
                f"{np.sqrt(np.mean(coefficient_sigma**2)):.2g}/"
                f"{np.max(coefficient_sigma):.2g}",
            )
        ),
        va="top",
        linespacing=1.8,
    )
    info_axis.legend(
        handles=(
            Patch(facecolor="0.75", edgecolor="0.55", label="CT18 90% / CT18 central"),
            Patch(
                facecolor="none",
                edgecolor="tab:orange",
                hatch="//////",
                label="projected 90% / CT18 central",
            ),
            Line2D([], [], color="black", label="CT18 central"),
            Line2D(
                [],
                [],
                color="tab:red",
                ls="--",
                label="projected central / CT18 central",
            ),
        ),
        loc="lower left",
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    ratio_bands_png = args.output_directory / "ct18nnlo_pod_ratio_bands.png"
    ratio_bands_pdf = args.output_directory / "ct18nnlo_pod_ratio_bands.pdf"
    fig.savefig(ratio_bands_png, dpi=180)
    fig.savefig(ratio_bands_pdf)
    plt.close(fig)

    for output in (
        band_png,
        band_pdf,
        ratio_png,
        ratio_pdf,
        ratio_bands_png,
        ratio_bands_pdf,
    ):
        print(f"Wrote {output.resolve()}")


if __name__ == "__main__":
    main()
