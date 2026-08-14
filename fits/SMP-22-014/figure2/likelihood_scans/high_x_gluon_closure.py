#!/usr/bin/env python3
"""Test central fitted-PDF closure of the full POD basis at high-x."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import lhapdf
import matplotlib.pyplot as plt
import numpy as np
import yaml

from evaluate_direct import parse_likelihood, render_parameters, render_steering
from scan_tools import read_lhagrid_first_q, run_xfitter, sha256, write_fixed_nuisances


def solve_coefficients(
    reference: np.ndarray, shifts: np.ndarray, target: np.ndarray
) -> np.ndarray:
    matrix = shifts.reshape(shifts.shape[0], -1).T
    displacement = (target - reference).reshape(-1)
    return np.linalg.solve(matrix.T @ matrix, matrix.T @ displacement)


def interval_metrics(
    x: np.ndarray, target: np.ndarray, reconstructed: np.ndarray
) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    scale = np.max(np.abs(target))
    for low, high in ((0.1, 0.3), (0.3, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 0.99)):
        mask = (x >= low) & (x < high) & (np.abs(target) > 1.0e-12 * scale)
        relative = reconstructed[mask] / target[mask] - 1.0
        result[f"{low:g}_to_{high:g}"] = {
            "point_count": int(mask.sum()),
            "rms_relative": float(np.sqrt(np.mean(relative * relative))),
            "maximum_absolute_relative": float(np.max(np.abs(relative))),
        }
    return result


def solve_high_x_weighted(
    reference: np.ndarray,
    shifts: np.ndarray,
    target: np.ndarray,
    x_grid: np.ndarray,
    gluon_index: int,
    relative_weight: float,
    relative_x_range: tuple[float, float],
) -> np.ndarray:
    """Augment dist0 with a relative-gluon closure term for a capability test."""
    absolute_matrix = shifts.reshape(shifts.shape[0], -1).T
    absolute_displacement = (target - reference).reshape(-1)
    gluon_target = target[gluon_index]
    low, high = relative_x_range
    mask = (
        (x_grid >= low)
        & (x_grid <= high)
        & (np.abs(gluon_target) > 1.0e-14 * np.max(np.abs(gluon_target)))
    )
    relative_matrix = (
        shifts[:, gluon_index, mask].T / gluon_target[mask, np.newaxis]
    )
    relative_displacement = (
        gluon_target[mask] - reference[gluon_index, mask]
    ) / gluon_target[mask]
    matrix = np.vstack(
        [absolute_matrix, relative_weight * relative_matrix]
    )
    displacement = np.concatenate(
        [absolute_displacement, relative_weight * relative_displacement]
    )
    return np.linalg.lstsq(matrix, displacement, rcond=None)[0]


def prepare_full_x_target(
    run_dir: Path,
    project_root: Path,
    figure2_dir: Path,
    scan_dir: Path,
    x_grid: np.ndarray,
    q_ext: float,
) -> Path:
    target_name = "direct_full_x_target"
    target_path = run_dir / "output" / target_name / f"{target_name}_0000.dat"
    if target_path.is_file():
        return target_path
    if run_dir.exists():
        raise RuntimeError(
            f"Incomplete target directory already exists: {run_dir}; move it aside before retrying"
        )

    reference_dir = figure2_dir / "reference_fit"
    fit = yaml.safe_load((reference_dir / "fit_result.yaml").read_text())
    nuisance_payload = yaml.safe_load((reference_dir / "nuisances.yaml").read_text())
    values = {
        parameter["name"]: float(parameter["value"])
        for parameter in fit["free_parameters"]
    }
    parameter_template = (
        reference_dir / "cards" / "covariance" / "parameters.yaml"
    ).read_text()
    parameters = render_parameters(
        parameter_template,
        values,
        {
            "name": target_name,
            "x_values": [float(value) for value in x_grid],
            "q_values": [
                q_ext,
                q_ext * (1.0 + 1.0e-6),
                q_ext * (1.0 + 2.0e-6),
            ],
        },
    )
    steering_template = (
        reference_dir / "cards" / "covariance" / "steering.txt"
    ).read_text()
    steering, _, _ = render_steering(steering_template)

    run_dir.mkdir(parents=True)
    (run_dir / "parameters.yaml").write_text(parameters)
    (run_dir / "steering.txt").write_text(steering)
    shutil.copy2(
        reference_dir / "cards" / "covariance" / "constants.yaml",
        run_dir / "constants.yaml",
    )
    write_fixed_nuisances(
        run_dir / "fixed_nuisances.dat",
        nuisance_payload["nuisances"],
        fit["fit_id"],
        nuisance_payload["stored_precision"],
    )
    (run_dir / "datafiles").symlink_to(
        project_root / "xfitter-datafiles", target_is_directory=True
    )
    (run_dir / "unpolarised.wgt").symlink_to(
        figure2_dir / "central" / "unpolarised.wgt"
    )
    run_xfitter(run_dir, project_root)
    likelihood = parse_likelihood(run_dir / "output" / "likelihood.txt")
    expected = float(fit["minimum"]["chi2"])
    closure = float(likelihood["total_chi2"]) - expected
    if abs(closure) > 0.01:
        raise RuntimeError(f"Direct full-x target failed likelihood closure: {closure:+.6g}")
    metadata = {
        "schema_version": 1,
        "fit_id": fit["fit_id"],
        "parameterization": "direct_HERAPDF",
        "q_ext_GeV": q_ext,
        "x_point_count": int(x_grid.size),
        "x_range": [float(x_grid[0]), float(x_grid[-1])],
        "total_chi2": float(likelihood["total_chi2"]),
        "reference_chi2": expected,
        "delta_chi2": closure,
        "member_file": str(target_path.relative_to(run_dir)),
    }
    (run_dir / "target_metadata.yaml").write_text(
        yaml.safe_dump(metadata, sort_keys=False, width=100)
    )
    return target_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("high_x_gluon_closure"),
        help="Output directory relative to this script (default: high_x_gluon_closure)",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    scan_dir = script_path.parent
    figure2_dir = scan_dir.parent
    project_root = script_path.parents[4]
    output = args.output if args.output.is_absolute() else scan_dir / args.output
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(project_root))
    from pod_projection.pod_projection import (
        LHAPDF_XGRID,
        QCD5_FLAVORS,
        ProjectionOperator,
    )

    lhapdf.setVerbosity(0)
    q_ext = 1.65
    full_x = LHAPDF_XGRID[36:]
    current_count = len(LHAPDF_XGRID[36:-20])
    target_path = prepare_full_x_target(
        output / "reference_direct_full_x",
        project_root,
        figure2_dir,
        scan_dir,
        full_x,
        q_ext,
    )
    target = read_lhagrid_first_q(
        target_path, QCD5_FLAVORS, full_x, q_ext
    )

    operator = ProjectionOperator.build(
        "250503_pod_basis_40k",
        100,
        QCD5_FLAVORS,
        full_x,
        q_ext,
        "dist0",
    )
    current_coefficients = solve_coefficients(
        operator.reference_grid[:, :current_count],
        operator.shift_grid[:, :, :current_count],
        target[:, :current_count],
    )
    full_coefficients = solve_coefficients(
        operator.reference_grid,
        operator.shift_grid,
        target,
    )
    current_reconstruction = operator.reference_grid + np.einsum(
        "i,ifx->fx", current_coefficients, operator.shift_grid
    )
    full_reconstruction = operator.reference_grid + np.einsum(
        "i,ifx->fx", full_coefficients, operator.shift_grid
    )
    gluon_index = QCD5_FLAVORS.index(21)
    relative_weight = 0.1
    relative_x_range = (0.05, 0.99)
    high_x_weighted_coefficients = solve_high_x_weighted(
        operator.reference_grid,
        operator.shift_grid,
        target,
        full_x,
        gluon_index,
        relative_weight,
        relative_x_range,
    )
    high_x_weighted_reconstruction = operator.reference_grid + np.einsum(
        "i,ifx->fx", high_x_weighted_coefficients, operator.shift_grid
    )

    stored_path = (
        scan_dir
        / "evaluations"
        / "reference_direct_with_projection"
        / "full_pod_projection.npz"
    )
    with np.load(stored_path, allow_pickle=False) as stored:
        stored_coefficients = np.asarray(stored["coefficients"], dtype=float)
    coefficient_reproduction = float(
        np.linalg.norm(current_coefficients - stored_coefficients)
        / np.linalg.norm(stored_coefficients)
    )

    target_gluon = target[gluon_index]
    current_gluon = current_reconstruction[gluon_index]
    full_gluon = full_reconstruction[gluon_index]
    weighted_gluon = high_x_weighted_reconstruction[gluon_index]
    cutoff = float(full_x[current_count - 1])
    metrics = {
        "current_projection_range": interval_metrics(
            full_x, target_gluon, current_gluon
        ),
        "full_x_projection_range": interval_metrics(
            full_x, target_gluon, full_gluon
        ),
        "high_x_weighted_projection": interval_metrics(
            full_x, target_gluon, weighted_gluon
        ),
    }

    np.savez_compressed(
        output / "high_x_gluon_closure.npz",
        x_grid=full_x,
        flavors=np.asarray(QCD5_FLAVORS),
        target_grid=target,
        pod_reference_grid=operator.reference_grid,
        pod_shift_grid=operator.shift_grid,
        current_coefficients=current_coefficients,
        full_x_coefficients=full_coefficients,
        high_x_weighted_coefficients=high_x_weighted_coefficients,
        current_reconstruction_grid=current_reconstruction,
        full_x_reconstruction_grid=full_reconstruction,
        high_x_weighted_reconstruction_grid=high_x_weighted_reconstruction,
        current_residual_grid=target - current_reconstruction,
        full_x_residual_grid=target - full_reconstruction,
        high_x_weighted_residual_grid=target - high_x_weighted_reconstruction,
        q_ext_GeV=q_ext,
        current_projection_xmax=cutoff,
    )

    plot_mask = (full_x >= 0.05) & (full_x <= relative_x_range[1])
    ratio_mask = plot_mask & (np.abs(target_gluon) > 1.0e-12 * np.max(np.abs(target_gluon)))
    figure, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(8.2, 7.6),
        sharex=True,
        gridspec_kw={"height_ratios": [1.65, 1.0]},
        constrained_layout=True,
    )
    top.plot(full_x[plot_mask], target_gluon[plot_mask], color="black", lw=2.2, label="direct HERAPDF")
    top.plot(
        full_x[plot_mask],
        current_gluon[plot_mask],
        color="#d95f02",
        lw=1.8,
        ls="--",
        label=rf"POD projected to $x\leq {cutoff:.3f}$",
    )
    top.plot(
        full_x[plot_mask],
        full_gluon[plot_mask],
        color="#1b9e77",
        lw=1.8,
        ls=":",
        label=r"POD projected to $x=1$",
    )
    top.plot(
        full_x[plot_mask],
        weighted_gluon[plot_mask],
        color="#7570b3",
        lw=1.8,
        ls=(0, (5, 2, 1, 2)),
        label=r"POD with relative high-$x$ gluon term",
    )
    if (
        np.all(target_gluon[plot_mask] > 0)
        and np.all(current_gluon[plot_mask] > 0)
        and np.all(full_gluon[plot_mask] > 0)
        and np.all(weighted_gluon[plot_mask] > 0)
    ):
        top.set_yscale("log")
    else:
        positive = np.abs(target_gluon[plot_mask])
        top.set_yscale("symlog", linthresh=max(np.min(positive[positive > 0]), 1.0e-14))
    top.set_ylabel(r"$xg(x,Q=1.65\,\mathrm{GeV})$")
    top.legend(frameon=False)
    top.grid(alpha=0.2)
    for values, color, style, label in (
        (current_gluon, "#d95f02", "--", "current range"),
        (full_gluon, "#1b9e77", ":", "full-x range"),
        (weighted_gluon, "#7570b3", "-.", "relative high-x term"),
    ):
        bottom.plot(
            full_x[ratio_mask],
            100.0 * np.abs(values[ratio_mask] / target_gluon[ratio_mask] - 1.0),
            color=color,
            lw=1.8,
            ls=style,
            label=label,
        )
    for axis in (top, bottom):
        axis.axvline(cutoff, color="0.45", lw=1.0, ls="-.")
    bottom.set_yscale("log")
    bottom.set_xlim(0.05, relative_x_range[1])
    bottom.set_xlabel(r"$x$")
    bottom.set_ylabel("|POD / direct - 1| [%]")
    bottom.grid(alpha=0.2)
    bottom.legend(frameon=False)
    figure.suptitle("Full 100-mode POD closure of the fitted high-x gluon")
    figure.savefig(output / "high_x_gluon_closure.png", dpi=190)
    figure.savefig(output / "high_x_gluon_closure.pdf")
    plt.close(figure)

    metadata = {
        "schema_version": 1,
        "status": "complete",
        "basis_set": "250503_pod_basis_40k",
        "basis_members": [1, 100],
        "target": "central fitted direct HERAPDF",
        "uncertainties_projected": False,
        "q_ext_GeV": q_ext,
        "flavors": list(QCD5_FLAVORS),
        "metric": "dist0",
        "current_projection": {
            "x_slice": [36, -20],
            "x_point_count": current_count,
            "x_max": cutoff,
        },
        "full_x_projection": {
            "x_slice": [36, None],
            "x_point_count": int(full_x.size),
            "x_max": float(full_x[-1]),
        },
        "high_x_weighted_capability_test": {
            "base_metric": "dist0 on all full-x flavor nodes",
            "additional_term": "relative gluon residual",
            "relative_weight": relative_weight,
            "relative_x_range": list(relative_x_range),
            "purpose": "basis capability diagnostic; not yet the configured likelihood projection",
        },
        "current_coefficient_reproduction_relative": coefficient_reproduction,
        "coefficient_change_full_x_relative": float(
            np.linalg.norm(full_coefficients - current_coefficients)
            / np.linalg.norm(current_coefficients)
        ),
        "coefficient_l2_norms": {
            "current": float(np.linalg.norm(current_coefficients)),
            "full_x": float(np.linalg.norm(full_coefficients)),
            "high_x_weighted": float(np.linalg.norm(high_x_weighted_coefficients)),
        },
        "global_residual_over_target_shift": {
            "current": float(
                np.linalg.norm(target - current_reconstruction)
                / np.linalg.norm(target - operator.reference_grid)
            ),
            "full_x": float(
                np.linalg.norm(target - full_reconstruction)
                / np.linalg.norm(target - operator.reference_grid)
            ),
            "high_x_weighted": float(
                np.linalg.norm(target - high_x_weighted_reconstruction)
                / np.linalg.norm(target - operator.reference_grid)
            ),
        },
        "gluon_closure": metrics,
        "files": {
            "arrays": "high_x_gluon_closure.npz",
            "arrays_sha256": sha256(output / "high_x_gluon_closure.npz"),
            "plot_png": "high_x_gluon_closure.png",
            "plot_pdf": "high_x_gluon_closure.pdf",
            "target_grid": str(target_path.relative_to(output)),
            "script": str(script_path.relative_to(project_root)),
            "script_sha256": sha256(script_path),
        },
    }
    (output / "closure_metadata.yaml").write_text(
        yaml.safe_dump(metadata, sort_keys=False, width=105)
    )
    print(yaml.safe_dump(metadata, sort_keys=False, width=105))


if __name__ == "__main__":
    main()
