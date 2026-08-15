#!/usr/bin/env python3
"""Run resumable direct-versus-full-POD one-dimensional likelihood scans."""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
import sys
from pathlib import Path

import lhapdf
import numpy as np
import yaml

from scan_tools import read_lhagrid_first_q, sha256


def coordinate_tag(value: float) -> str:
    return f"sigma_{value:+.6f}".replace("+", "p").replace("-", "m").replace(".", "d")


def load_evaluation(path: Path) -> dict[str, object]:
    evaluation = path / "evaluation.yaml"
    if not evaluation.is_file():
        raise RuntimeError(f"Missing completed evaluation {evaluation}")
    result = yaml.safe_load(evaluation.read_text())
    if result.get("status") not in ("complete",):
        raise RuntimeError(f"Evaluation is not complete: {evaluation}")
    return result


def remove_incomplete_evaluation(path: Path) -> None:
    """Remove only a failed evaluator directory so a scan can resume safely."""
    if path.exists() and not (path / "evaluation.yaml").is_file():
        shutil.rmtree(path)


def group_totals(evaluation: dict[str, object]) -> dict[str, float]:
    groups = evaluation["likelihood_groups"]
    hera = float(groups["HERA"]["total_chi2"])
    cms = float(groups["CMS"]["total_chi2"])
    return {"HERA": hera, "CMS": cms, "sum": hera + cms}


def run_logged(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        process = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if process.returncode:
        tail = "\n".join(log_path.read_text().splitlines()[-30:])
        raise RuntimeError(
            f"Command failed with exit code {process.returncode}: {' '.join(command)}\n{tail}"
        )


def save_results(
    output: Path,
    rows: list[dict[str, object]],
    fit: dict[str, object],
    covariance: dict[str, np.ndarray],
    config: dict[str, object],
    operator: object,
    started: dt.datetime,
    scan_dir: Path,
) -> None:
    ordered_names = [str(name) for name in covariance["parameter_names"]]
    arrays: dict[str, np.ndarray] = {
        "parameter_name": np.asarray([row["parameter_name"] for row in rows]),
        "parameter_index": np.asarray([row["parameter_index"] for row in rows], dtype=int),
        "scan_coordinate_sigma": np.asarray([row["coordinate"] for row in rows], dtype=float),
        "parameter_value": np.asarray([row["parameter_value"] for row in rows], dtype=float),
        "evaluation_status": np.asarray(["complete" for row in rows]),
        "direct_chi2_HERA": np.asarray([row["direct"]["HERA"] for row in rows]),
        "direct_chi2_CMS": np.asarray([row["direct"]["CMS"] for row in rows]),
        "direct_chi2_sum": np.asarray([row["direct"]["sum"] for row in rows]),
        "full_pod_chi2_HERA": np.asarray([row["pod"]["HERA"] for row in rows]),
        "full_pod_chi2_CMS": np.asarray([row["pod"]["CMS"] for row in rows]),
        "full_pod_chi2_sum": np.asarray([row["pod"]["sum"] for row in rows]),
        "delta_chi2_HERA": np.asarray([row["pod"]["HERA"] - row["direct"]["HERA"] for row in rows]),
        "delta_chi2_CMS": np.asarray([row["pod"]["CMS"] - row["direct"]["CMS"] for row in rows]),
        "delta_chi2_sum": np.asarray([row["pod"]["sum"] - row["direct"]["sum"] for row in rows]),
        "projection_relative_residual": np.asarray([row["projection_relative_residual"] for row in rows]),
        "external_parameter_values": np.asarray([row["external_values"] for row in rows], dtype=float),
        "full_pod_coefficients": np.asarray([row["coefficients"] for row in rows], dtype=float),
        "direct_evaluation_directory": np.asarray([row["direct_directory"] for row in rows]),
        "pod_evaluation_directory": np.asarray([row["pod_directory"] for row in rows]),
    }
    np.savez_compressed(output / "scan_results.npz", **arrays)
    metadata = {
        "schema_version": 1,
        "status": "partial" if len(rows) == 0 else "complete_points_stored",
        "reference_fit_id": fit["fit_id"],
        "created_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "started_utc": started.replace(microsecond=0).isoformat(),
        "row_count": len(rows),
        "parameter_order": ordered_names,
        "minimization": False,
        "profiling": False,
        "nuisance_treatment": "fixed at stored global-best-fit values",
        "projection": {
            **config["projection"],
            "matrix_rank": operator.rank,
            "normal_matrix_condition_number": operator.condition_number,
        },
        "files": {
            "results": "scan_results.npz",
            "results_sha256": sha256(output / "scan_results.npz"),
            "scan_config": str((scan_dir / "scan_config.yaml").resolve()),
            "scan_config_sha256": sha256(scan_dir / "scan_config.yaml"),
        },
    }
    (output / "scan_metadata.yaml").write_text(
        yaml.safe_dump(metadata, sort_keys=False, width=110)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--parameters",
        help="Comma-separated parameter names (default: all 16 independent parameters)",
    )
    parser.add_argument(
        "--coordinates",
        help="Comma-separated HESSE-sigma coordinates (default: scan_config range/points)",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    scan_dir = script_path.parent
    figure2_dir = scan_dir.parent
    project_root = script_path.parents[4]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load((scan_dir / "scan_config.yaml").read_text())
    fit = yaml.safe_load((figure2_dir / "reference_fit" / "fit_result.yaml").read_text())
    with np.load(figure2_dir / "reference_fit" / "covariance.npz", allow_pickle=False) as payload:
        covariance = {key: np.array(payload[key]) for key in payload.files}
    names = [str(name) for name in covariance["parameter_names"]]
    selected = (
        [name.strip() for name in args.parameters.split(",") if name.strip()]
        if args.parameters
        else names
    )
    unknown = sorted(set(selected) - set(names))
    if unknown:
        raise SystemExit(f"Unknown scan parameters: {unknown}")
    if args.coordinates:
        coordinates = [float(value) for value in args.coordinates.split(",")]
    else:
        low, high = config["scan"]["range"]
        coordinates = np.linspace(low, high, int(config["scan"]["points"])).tolist()
    if len(set(coordinates)) != len(coordinates):
        raise SystemExit("Duplicate scan coordinates are not allowed")

    lhapdf.setVerbosity(0)
    sys.path.insert(0, str(project_root))
    from pod_projection.pod_projection import LHAPDF_XGRID
    from projection_metrics import Figure2ProjectionOperator

    projection = config["projection"]
    start, stop = projection["x_slice"]
    x_grid = LHAPDF_XGRID[start:stop]
    flavors = tuple(int(pid) for pid in projection["flavors"])
    q_ext = float(projection["q_ext_GeV"])
    operator = Figure2ProjectionOperator.build(
        projection["basis_set"],
        int(projection["coefficient_count"]),
        flavors,
        x_grid,
        q_ext,
        projection["metric"],
        relative_weight=float(projection.get("relative_weight", 0.1)),
        relative_x_range=tuple(projection.get("relative_x_range", [0.05, 0.99])),
        relative_floor=float(projection.get("relative_floor", 1.0e-12)),
        relative_valence_weight=float(projection.get("relative_valence_weight", 0.0)),
        relative_valence_x_range=tuple(
            projection.get("relative_valence_x_range", [1.0e-4, 0.1])
        ),
        relative_valence_floor=float(
            projection.get("relative_valence_floor", 1.0e-12)
        ),
        relative_f2_weight=float(projection.get("relative_f2_weight", 0.0)),
        relative_f2_x_range=tuple(
            projection.get("relative_f2_x_range", [1.0e-4, 0.1])
        ),
        relative_f2_floor=float(projection.get("relative_f2_floor", 1.0e-12)),
    )

    values = np.asarray(covariance["values"], dtype=float)
    errors = np.asarray(covariance["hesse_errors"], dtype=float)
    reference_alpha = float(values[names.index("alphas")])
    # Keep one central evaluation inside each scan output. This makes it
    # reusable by all parameter axes while tying it to the exact projection
    # contract used for that scan.
    canonical_direct = output / "runs" / "_reference" / "direct"
    canonical_pod = output / "runs" / "_reference" / "full_pod"
    rows: list[dict[str, object]] = []
    started = dt.datetime.now(dt.timezone.utc)

    for parameter in selected:
        index = names.index(parameter)
        for coordinate in sorted(coordinates):
            point_value = float(values[index] + coordinate * errors[index])
            point_root = output / "runs" / parameter / coordinate_tag(coordinate)
            if coordinate == 0.0:
                work_root = canonical_direct.parent
                direct_dir = canonical_direct
                pod_dir = canonical_pod
                projection_path = work_root / "full_pod_projection.npz"
            else:
                work_root = point_root
                direct_dir = point_root / "direct"
                pod_dir = point_root / "full_pod"
                projection_path = point_root / "full_pod_projection.npz"
            print(
                f"[{len(rows) + 1}/{len(selected) * len(coordinates)}] "
                f"{parameter} {coordinate:+.3f} sigma -> {point_value:.9g}",
                flush=True,
            )

            direct_source = direct_dir
            pod_source = pod_dir
            projection_source = projection_path
            if not (direct_dir / "evaluation.yaml").is_file():
                remove_incomplete_evaluation(direct_dir)
                direct_arguments = [
                    sys.executable,
                    str(scan_dir / "evaluate_direct.py"),
                ]
                if coordinate != 0.0:
                    direct_arguments.extend(
                        ["--parameter", f"{parameter}={point_value:.17g}"]
                    )
                direct_arguments.extend(
                    ["--export-projection-pdf", "--output", str(direct_dir)]
                )
                run_logged(
                    direct_arguments,
                    work_root / "evaluate_direct.console.log",
                )
            direct_evaluation = load_evaluation(direct_dir)
            target_path = direct_dir / str(
                direct_evaluation["projection_pdf"]["member_file"]
            )
            target = read_lhagrid_first_q(target_path, flavors, x_grid, q_ext)
            projected, coefficients, residual = operator.project_grid(target)
            denominator = np.linalg.norm(target - operator.reference_grid)
            relative_residual = float(
                np.linalg.norm(residual) / denominator if denominator else 0.0
            )
            work_root.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                projection_path,
                coefficients=coefficients,
                target_grid=target,
                projected_grid=projected,
                residual_grid=residual,
                x_grid=x_grid,
                flavors=np.asarray(flavors),
                q_ext_GeV=q_ext,
                projection_relative_residual=relative_residual,
            )
            if not (pod_dir / "evaluation.yaml").is_file():
                remove_incomplete_evaluation(pod_dir)
                alpha = point_value if parameter == "alphas" else reference_alpha
                run_logged(
                    [
                        sys.executable,
                        str(scan_dir / "evaluate_pod.py"),
                        "--coefficients",
                        str(projection_path),
                        "--alphas",
                        f"{alpha:.17g}",
                        "--source-direct-evaluation",
                        str(direct_dir),
                        "--output",
                        str(pod_dir),
                    ],
                    work_root / "evaluate_pod.console.log",
                )

            direct_evaluation = load_evaluation(direct_source)
            pod_evaluation = load_evaluation(pod_source)
            with np.load(projection_source, allow_pickle=False) as payload:
                coefficients = np.asarray(payload["coefficients"], dtype=float)
                target = np.asarray(payload["target_grid"], dtype=float)
                residual = np.asarray(payload["residual_grid"], dtype=float)
                if "projection_relative_residual" in payload:
                    relative_residual = float(payload["projection_relative_residual"])
                else:
                    denominator = np.linalg.norm(target - operator.reference_grid)
                    relative_residual = float(np.linalg.norm(residual) / denominator)
            external_values = values.copy()
            external_values[index] = point_value
            rows.append(
                {
                    "parameter_name": parameter,
                    "parameter_index": index,
                    "coordinate": coordinate,
                    "parameter_value": point_value,
                    "external_values": external_values,
                    "coefficients": coefficients,
                    "projection_relative_residual": relative_residual,
                    "direct": group_totals(direct_evaluation),
                    "pod": group_totals(pod_evaluation),
                    "direct_directory": str(direct_source.resolve()),
                    "pod_directory": str(pod_source.resolve()),
                }
            )
            save_results(
                output, rows, fit, covariance, config, operator, started, scan_dir
            )

    metadata_path = output / "scan_metadata.yaml"
    metadata = yaml.safe_load(metadata_path.read_text())
    metadata["status"] = "complete"
    metadata["completed_utc"] = dt.datetime.now(dt.timezone.utc).replace(
        microsecond=0
    ).isoformat()
    metadata["requested_parameters"] = selected
    metadata["requested_coordinates"] = sorted(coordinates)
    metadata_path.write_text(yaml.safe_dump(metadata, sort_keys=False, width=110))
    print(f"Stored {len(rows)} points in {output / 'scan_results.npz'}")


if __name__ == "__main__":
    main()
