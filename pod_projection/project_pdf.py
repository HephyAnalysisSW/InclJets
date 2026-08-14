#!/usr/bin/env python3
"""Project every member of an external LHAPDF set onto the native POD basis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lhapdf
import numpy as np

try:
    from .native_pod_basis_40k import PDF_SET, pdfset_size
    from .pod_projection import (
        LHAPDF_XGRID,
        METRICS,
        Q0,
        ProjectionOperator,
        coefficient_covariance,
        hessian_symmetric_band,
        parse_flavors,
    )
except ImportError:  # Support `python project_pdf.py` from this directory.
    from native_pod_basis_40k import PDF_SET, pdfset_size
    from pod_projection import (
        LHAPDF_XGRID,
        METRICS,
        Q0,
        ProjectionOperator,
        coefficient_covariance,
        hessian_symmetric_band,
        parse_flavors,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project an external LHAPDF set onto native POD shifts."
    )
    parser.add_argument("--basis-set", default=PDF_SET)
    parser.add_argument("--target-set", default="CT18NNLO")
    parser.add_argument("--q", type=float, default=Q0)
    parser.add_argument("--flavors", default="qcd5")
    parser.add_argument("--n-basis", type=int, default=100)
    parser.add_argument("--x-start", type=int, default=36)
    parser.add_argument("--x-stop", type=int, default=-20)
    parser.add_argument("--metric", choices=METRICS, default="dist0")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent / "outputs" / "ct18nnlo_pod_projection.npz",
    )
    return parser.parse_args()


def weighted_norm(values: np.ndarray, weights: np.ndarray) -> float:
    flat = np.asarray(values, dtype=float).reshape(-1)
    return float(np.sqrt(max(float(flat @ (weights * flat)), 0.0)))


def main() -> None:
    args = parse_args()
    try:
        lhapdf.setVerbosity(0)
    except AttributeError:
        pass

    target_info = lhapdf.getPDFSet(args.target_set)
    n_members = pdfset_size(args.target_set)
    error_type_value = target_info.errorType
    error_type = str(
        error_type_value() if callable(error_type_value) else error_type_value
    )
    try:
        confidence_level = target_info.get_entry("ErrorConfLevel")
        confidence_level = float(confidence_level)
        if confidence_level.is_integer():
            confidence_level = int(confidence_level)
    except Exception:
        confidence_level = None
    x_grid = np.asarray(LHAPDF_XGRID[args.x_start : args.x_stop], dtype=float)
    if not len(x_grid):
        raise ValueError("The requested x-grid slice is empty")
    flavors = parse_flavors(args.flavors, args.basis_set, args.target_set)

    print(
        f"Building {args.n_basis}-mode operator on "
        f"{len(flavors)} flavors x {len(x_grid)} x-points ..."
    )
    operator = ProjectionOperator.build(
        basis_set=args.basis_set,
        n_basis=args.n_basis,
        flavors=flavors,
        x_grid=x_grid,
        q=args.q,
        metric=args.metric,
    )

    target_grids = []
    projected_grids = []
    residual_grids = []
    coefficients = []
    relative_residuals = []
    for member in range(n_members):
        target = operator.evaluate(lhapdf.mkPDF(args.target_set, member))
        projected, coeffs, residual = operator.project_grid(target)
        target_grids.append(target)
        projected_grids.append(projected)
        residual_grids.append(residual)
        coefficients.append(coeffs)
        denominator = weighted_norm(target - operator.reference_grid, operator.weights)
        numerator = weighted_norm(residual, operator.weights)
        relative_residuals.append(numerator / denominator if denominator else np.nan)
        if member == 0 or (member + 1) % 10 == 0 or member + 1 == n_members:
            print(f"  projected member {member:02d}/{n_members - 1:02d}")

    target_grids = np.asarray(target_grids)
    projected_grids = np.asarray(projected_grids)
    residual_grids = np.asarray(residual_grids)
    coefficients = np.asarray(coefficients)
    relative_residuals = np.asarray(relative_residuals)

    target_band = np.empty((0,), dtype=float)
    projected_band = np.empty((0,), dtype=float)
    coefficient_cov = np.empty((0, 0), dtype=float)
    coefficient_displacements = np.empty((0, args.n_basis), dtype=float)
    hessian_pairs = 0
    if (
        "hessian" in error_type.lower()
        and n_members >= 3
        and (n_members - 1) % 2 == 0
    ):
        target_band = hessian_symmetric_band(target_grids)
        projected_band = hessian_symmetric_band(projected_grids)
        coefficient_cov, coefficient_displacements = coefficient_covariance(coefficients)
        hessian_pairs = (n_members - 1) // 2

    covariance_rank = (
        int(np.linalg.matrix_rank(coefficient_cov)) if coefficient_cov.size else 0
    )
    metadata = {
        "basis_set": args.basis_set,
        "target_set": args.target_set,
        "q_GeV": args.q,
        "metric": args.metric,
        "n_basis": args.n_basis,
        "n_members": n_members,
        "target_error_type": error_type,
        "x_start": args.x_start,
        "x_stop": args.x_stop,
        "hessian_pairs": hessian_pairs,
        "hessian_confidence_level": confidence_level,
        "convention": "LHAPDF x*f(x,Q)",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
        x_grid=x_grid,
        flavors=np.asarray(flavors, dtype=int),
        reference_grid=operator.reference_grid,
        basis_shifts=operator.shift_grid,
        target_grids=target_grids,
        projected_grids=projected_grids,
        residual_grids=residual_grids,
        coefficients=coefficients,
        relative_residuals=relative_residuals,
        target_hessian_band=target_band,
        projected_hessian_band=projected_band,
        coefficient_covariance=coefficient_cov,
        coefficient_displacements=coefficient_displacements,
        gram=operator.gram,
    )

    finite_residuals = relative_residuals[np.isfinite(relative_residuals)]
    print(f"Rows / modes              : {operator.matrix.shape[0]} / {args.n_basis}")
    print(f"rank(X)                   : {operator.rank}")
    print(f"cond(X^T W X)             : {operator.condition_number:.6e}")
    print(f"central residual / shift  : {relative_residuals[0]:.6e}")
    print(
        "all-member residual ratio: "
        f"min/median/max = {np.min(finite_residuals):.3e} / "
        f"{np.median(finite_residuals):.3e} / {np.max(finite_residuals):.3e}"
    )
    if hessian_pairs:
        print(f"Hessian pairs             : {hessian_pairs}")
        print(f"rank(coefficient cov.)    : {covariance_rank}/{args.n_basis}")
    print(f"Wrote                     : {args.output.resolve()}")


if __name__ == "__main__":
    main()
