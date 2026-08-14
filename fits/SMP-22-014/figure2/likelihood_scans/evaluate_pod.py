#!/usr/bin/env python3
"""Evaluate the fixed-nuisance likelihood for one full-POD coefficient vector."""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
from pathlib import Path

import numpy as np
import yaml

from evaluate_direct import render_steering
from scan_tools import (
    likelihood_groups,
    parse_likelihood,
    run_xfitter,
    sha256,
    xfitter_library,
    write_checksums,
    write_fixed_nuisances,
)


def load_coefficients(path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as payload:
        for key in ("coefficients", "full_pod_coefficients"):
            if key in payload:
                coefficients = np.asarray(payload[key], dtype=float).reshape(-1)
                break
        else:
            raise ValueError(f"{path} has no coefficients array")
    if coefficients.shape != (100,) or not np.all(np.isfinite(coefficients)):
        raise ValueError("Expected exactly 100 finite full-POD coefficients")
    return coefficients


def render_parameters(coefficients: np.ndarray, alphas: float) -> str:
    lines = [
        "# Full 100-mode POD point evaluation; all parameters are fixed.",
        "Minimizer: MINUIT",
        "MINUIT:",
        "  Commands: |",
        "    call fcn 3",
        "  doErrors: None",
        "",
        "Parameters:",
    ]
    lines.extend(
        f"  pod_c_{index}: [{value:.17g}, 0.0]"
        for index, value in enumerate(coefficients, start=1)
    )
    lines.extend(
        [
            f"  alphas: [{alphas:.17g}, 0.0]",
            "",
            "DefaultDecomposition: proton",
            "Decompositions:",
            "  proton:",
            "    ? !include pod.yaml",
            "",
            "DefaultEvolution: proton-QCDNUM",
            "Evolutions:",
            "  proton-QCDNUM:",
            "    ? !include evolutions/QCDNUM.yaml",
            "    decomposition: proton",
            "",
            "Order: NNLO",
            "NFlavour: 5",
            "isFFNS: 0",
            "Q0: 1.65",
            "",
            "? !include constants.yaml",
            "",
            "byReaction:",
            "  RT_DISNC:",
            "    ? !include reactions/RT_DISNC.yaml",
            "",
            "hf_scheme_DISNC:",
            "  defaultValue: RT_DISNC",
            "",
            "hf_scheme_DISCC:",
            "  defaultValue: BaseDISCC",
            "",
            "OutputDirectory: output",
        ]
    )
    return "\n".join(lines) + "\n"


def render_pod_card() -> str:
    return """# Immutable full-basis POD definition for projection-bias scans.
class: POD
set: 250503_pod_basis_40k
reference_member: 0
basis_members: all
active_flavors: all
coefficient_prefix: pod_c_
q0: 1.65
require_matching_q0: true
value_convention: xfx
member_convention: member_minus_reference

projection:
  grid: LHAPDF_XGRID
  x_slice: [36, null]
  flavors: [21, 2, -2, 1, -1, 3, -3, 4, -4, 5, -5]
  # Absolute residual plus relative gluon residual; diagnostic projection
  # metric only, not an additional likelihood or prior.
  metric: relative_gluon
  relative_weight: 0.1
  relative_x_range: [0.05, 0.99]
  relative_floor: 1.0e-12
  # Relative u_v and d_v closure terms are used during external projection;
  # they are not extra xFitter likelihood terms or PDF priors.
  relative_valence_weight: 0.1
  relative_valence_x_range: [1.0e-4, 0.1]
  relative_valence_floor: 1.0e-12
  # Charge-weighted photon-exchange F2 proxy used during projection only.
  relative_f2_weight: 3.0
  relative_f2_x_range: [1.0e-4, 0.1]
  relative_f2_floor: 1.0e-12
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", required=True, type=Path)
    parser.add_argument("--alphas", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-direct-evaluation", type=Path)
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    scan_dir = script_path.parent
    figure2_dir = scan_dir.parent
    project_root = script_path.parents[4]
    reference_dir = figure2_dir / "reference_fit"
    run_dir = args.output.resolve()
    if run_dir.exists():
        raise SystemExit(f"Refusing to overwrite evaluation directory {run_dir}")

    coefficients = load_coefficients(args.coefficients)
    fit = yaml.safe_load((reference_dir / "fit_result.yaml").read_text())
    nuisance_payload = yaml.safe_load((reference_dir / "nuisances.yaml").read_text())
    nuisances = nuisance_payload["nuisances"]
    steering_template = (
        reference_dir / "cards" / "covariance" / "steering.txt"
    ).read_text()
    steering, input_files, corr_files = render_steering(steering_template)

    run_dir.mkdir(parents=True)
    shutil.copy2(script_path, run_dir / script_path.name)
    shutil.copy2(args.coefficients, run_dir / "coefficients.npz")
    (run_dir / "parameters.yaml").write_text(
        render_parameters(coefficients, args.alphas)
    )
    (run_dir / "pod.yaml").write_text(render_pod_card())
    (run_dir / "steering.txt").write_text(steering)
    shutil.copy2(
        reference_dir / "cards" / "covariance" / "constants.yaml",
        run_dir / "constants.yaml",
    )
    write_fixed_nuisances(
        run_dir / "fixed_nuisances.dat",
        nuisances,
        fit["fit_id"],
        nuisance_payload["stored_precision"],
    )
    (run_dir / "datafiles").symlink_to(
        project_root / "xfitter-datafiles", target_is_directory=True
    )
    (run_dir / "unpolarised.wgt").symlink_to(
        figure2_dir / "central" / "unpolarised.wgt"
    )

    started = dt.datetime.now(dt.timezone.utc)
    run_xfitter(run_dir, project_root)
    finished = dt.datetime.now(dt.timezone.utc)
    likelihood = parse_likelihood(run_dir / "output" / "likelihood.txt")
    if likelihood["free_parameter_count"] != 0:
        raise SystemExit("POD evaluation unexpectedly had free parameters")
    if likelihood["nuisance_treatment"] != "fixed":
        raise SystemExit("POD evaluation unexpectedly profiled nuisances")
    groups = likelihood_groups(likelihood, nuisances)

    binary = project_root / "install" / "xfitter" / "bin" / "xfitter"
    library = xfitter_library(project_root)
    source_files = [
        project_root / "xfitter" / "pdfdecomps" / "POD" / "PODPdfDecomposition.cc",
        project_root / "xfitter" / "pdfdecomps" / "POD" / "PODBasis.cc",
        project_root / "xfitter" / "include" / "steering.inc",
        project_root / "xfitter" / "src" / "GetChisquare.f",
        project_root / "xfitter" / "src" / "fcn.f",
    ]
    result = {
        "schema_version": 1,
        "status": "complete",
        "fit_id": fit["fit_id"],
        "parameterization": "full_100_mode_POD",
        "basis_set": "250503_pod_basis_40k",
        "basis_members": [1, 100],
        "coefficient_count": 100,
        "alphas": args.alphas,
        "minimization": False,
        "profiling": False,
        "nuisance_source": str(
            (reference_dir / "nuisances.yaml").relative_to(project_root)
        ),
        "nuisance_count": len(nuisances),
        "input_files": input_files,
        "correlation_files": corr_files,
        "source_direct_evaluation": (
            str(args.source_direct_evaluation.resolve())
            if args.source_direct_evaluation
            else None
        ),
        "coefficients": {
            "file": "coefficients.npz",
            "sha256": sha256(run_dir / "coefficients.npz"),
            "minimum": float(np.min(coefficients)),
            "maximum": float(np.max(coefficients)),
            "l2_norm": float(np.linalg.norm(coefficients)),
        },
        "likelihood": likelihood,
        "likelihood_groups": groups,
        "runtime": {
            "started_utc": started.replace(microsecond=0).isoformat(),
            "finished_utc": finished.replace(microsecond=0).isoformat(),
            "wall_seconds": (finished - started).total_seconds(),
        },
        "software": {
            "xfitter_binary_sha256": sha256(binary),
            "xfitter_library_sha256": sha256(library),
            "xfitter_commit": subprocess.check_output(
                ["git", "-C", str(project_root / "xfitter"), "rev-parse", "HEAD"],
                text=True,
            ).strip(),
            "source_sha256": {
                str(path.relative_to(project_root)): sha256(path)
                for path in source_files
            },
        },
    }
    (run_dir / "evaluation.yaml").write_text(
        yaml.safe_dump(result, sort_keys=False, width=100)
    )
    write_checksums(run_dir)
    print(f"Evaluation stored in {run_dir}")
    print(f"total chi2 = {float(likelihood['total_chi2']):.12f}")


if __name__ == "__main__":
    main()
