#!/usr/bin/env python3
"""Evaluate the direct xFitter likelihood with no minimization or profiling."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from scan_tools import xfitter_environment, xfitter_library


DATASET_GROUPS = {
    "HERA": range(1, 8),
    "CMS": range(8, 12),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_overrides(items: list[str], allowed: set[str]) -> dict[str, float]:
    result: dict[str, float] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Parameter override must be NAME=VALUE, got {item!r}")
        name, value = item.split("=", 1)
        if name not in allowed:
            raise ValueError(f"Unknown independent parameter {name!r}")
        if name in result:
            raise ValueError(f"Duplicate override for {name!r}")
        result[name] = float(value)
    return result


def render_parameters(
    template: str,
    values: dict[str, float],
    projection_export: dict[str, object] | None = None,
) -> str:
    minuit_start = template.index("MINUIT:")
    parameters_start = template.index("\nParameters:", minuit_start)
    template = (
        template[:minuit_start]
        + "MINUIT:\n  Commands: |\n    call fcn 3\n  doErrors: None\n"
        + template[parameters_start:]
    )

    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
    for name, value in values.items():
        pattern = re.compile(
            rf"(?m)^(  {re.escape(name)}:\s*\[\s*){number}(\s*,\s*){number}(\s*\])"
        )
        template, count = pattern.subn(
            rf"\g<1>{value:.12g}\g<2>0.0\g<3>", template, count=1
        )
        if count != 1:
            raise RuntimeError(f"Could not fix parameter {name!r} in template")
    if projection_export is not None:
        if "WriteLHAPDF6:" in template:
            raise RuntimeError("Template already contains a WriteLHAPDF6 block")
        export_block = {
            "WriteLHAPDF6": {
                "name": projection_export["name"],
                "evolution": "proton-QCDNUM",
                "description": "Direct HERAPDF scan point on the POD projection grid",
                "Xvalues": projection_export["x_values"],
                "Qvalues": projection_export["q_values"],
            }
        }
        template = template.rstrip() + "\n\n" + yaml.safe_dump(
            export_block, sort_keys=False, width=120
        )
    return template


def render_steering(
    template: str, fixed_nuisances: bool = True
) -> tuple[str, list[str], list[str]]:
    xfitter_match = re.search(r"&xFitter\b.*?&End", template, flags=re.DOTALL)
    if not xfitter_match:
        raise RuntimeError("Missing &xFitter block in steering template")
    if fixed_nuisances:
        xfitter_block = xfitter_match.group(0).replace(
            "&End",
            "  UseFixedNuisances = True\n"
            "  FixedNuisanceFile = 'fixed_nuisances.dat'\n"
            "&End",
        )
        template = (
            template[: xfitter_match.start()]
            + xfitter_block
            + template[xfitter_match.end() :]
        )

    input_match = re.search(r"&InFiles\b.*?&End", template, flags=re.DOTALL)
    corr_match = re.search(r"&InCorr\b.*?&End", template, flags=re.DOTALL)
    if not input_match or not corr_match:
        raise RuntimeError("Missing &InFiles or &InCorr block in steering template")
    all_inputs = re.findall(r"'([^']+)'", input_match.group(0))
    all_corr = re.findall(r"'([^']+)'", corr_match.group(0))
    if len(all_inputs) != 11 or len(all_corr) != 10:
        raise RuntimeError(
            f"Expected 11 input and 10 correlation files; found {len(all_inputs)} and {len(all_corr)}"
        )

    return template, all_inputs, all_corr


def parse_likelihood(path: Path) -> dict[str, object]:
    integer_fields = {
        "schema_version",
        "data_point_count",
        "free_parameter_count",
        "nuisance_count",
        "dataset_count",
    }
    float_fields = {
        "total_chi2",
        "data_chi2",
        "correlated_penalty_chi2",
        "log_penalty_chi2",
        "additional_penalty_chi2",
    }
    result: dict[str, object] = {}
    for line in path.read_text().splitlines():
        key, value = line.split(maxsplit=1)
        if key in integer_fields:
            result[key] = int(value)
        elif key in float_fields or re.fullmatch(
            r"dataset_\d+_(?:data|log_penalty)_chi2", key
        ):
            result[key] = float(value)
        else:
            result[key] = value
    missing = (integer_fields | float_fields | {"nuisance_treatment"}) - result.keys()
    if missing:
        raise RuntimeError(f"Missing likelihood fields: {sorted(missing)}")
    return result


def write_checksums(run_dir: Path) -> None:
    payloads = sorted(
        path
        for path in run_dir.rglob("*")
        if path.is_file() and not path.is_symlink() and path.name != "files.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(run_dir)}" for path in payloads]
    (run_dir / "files.sha256").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parameter",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="Override one independent reference parameter (repeatable)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="New evaluation directory (default: evaluations/reference_direct_combined)",
    )
    parser.add_argument(
        "--export-projection-pdf",
        action="store_true",
        help="Also export the direct PDF on the exact full-POD projection grid",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    scan_dir = script_path.parent
    figure2_dir = scan_dir.parent
    project_root = script_path.parents[4]
    reference_dir = figure2_dir / "reference_fit"
    run_dir = (
        args.output
        or scan_dir / "evaluations" / "reference_direct_combined"
    ).resolve()
    if run_dir.exists():
        raise SystemExit(f"Refusing to overwrite evaluation directory {run_dir}")

    fit = yaml.safe_load((reference_dir / "fit_result.yaml").read_text())
    nuisance_payload = yaml.safe_load((reference_dir / "nuisances.yaml").read_text())
    reference_values = {
        parameter["name"]: float(parameter["value"])
        for parameter in fit["free_parameters"]
    }
    overrides = parse_overrides(args.parameter, set(reference_values))
    values = reference_values | overrides

    parameter_template = (
        reference_dir / "cards" / "covariance" / "parameters.yaml"
    ).read_text()
    steering_template = (
        reference_dir / "cards" / "covariance" / "steering.txt"
    ).read_text()
    projection_export = None
    if args.export_projection_pdf:
        sys.path.insert(0, str(project_root))
        from pod_projection.pod_projection import LHAPDF_XGRID

        scan_config = yaml.safe_load((scan_dir / "scan_config.yaml").read_text())
        projection = scan_config["projection"]
        start, stop = projection["x_slice"]
        q_ext = float(projection["q_ext_GeV"])
        projection_export = {
            "name": "direct_projection_target",
            "x_values": [float(x) for x in LHAPDF_XGRID[start:stop]],
            # Keep three nearby Q nodes for a valid LHAPDF interpolation grid. The
            # first is exactly Q_ext;
            # projection code reads that node directly, without interpolation.
            "q_values": [
                q_ext,
                q_ext * (1.0 + 1.0e-6),
                q_ext * (1.0 + 2.0e-6),
            ],
        }
    parameters = render_parameters(parameter_template, values, projection_export)
    steering, input_files, corr_files = render_steering(steering_template)
    selected_nuisances = nuisance_payload["nuisances"]

    run_dir.mkdir(parents=True)
    shutil.copy2(script_path, run_dir / "evaluate_direct.py")
    (run_dir / "parameters.yaml").write_text(parameters)
    (run_dir / "steering.txt").write_text(steering)
    shutil.copy2(
        reference_dir / "cards" / "covariance" / "constants.yaml",
        run_dir / "constants.yaml",
    )
    nuisance_lines = [
        "# local_index source_name fixed_shift",
        f"# source_fit_id {fit['fit_id']}",
        f"# stored_precision {nuisance_payload['stored_precision']}",
    ]
    nuisance_lines.extend(
        f"{local_index:4d} {nuisance['name']} {float(nuisance['shift']):.17g}"
        for local_index, nuisance in enumerate(selected_nuisances, start=1)
    )
    (run_dir / "fixed_nuisances.dat").write_text("\n".join(nuisance_lines) + "\n")
    (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
    (run_dir / "unpolarised.wgt").symlink_to(
        figure2_dir / "central" / "unpolarised.wgt"
    )

    xfitter_binary = project_root / "install" / "xfitter" / "bin" / "xfitter"
    xfitter_library_path = xfitter_library(project_root)

    started = dt.datetime.now(dt.timezone.utc)
    with (run_dir / "run.log").open("w") as log:
        process = subprocess.Popen(
            [str(xfitter_binary)],
            cwd=run_dir,
            env=xfitter_environment(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        return_code = process.wait()
    finished = dt.datetime.now(dt.timezone.utc)
    if return_code != 0:
        raise SystemExit(f"xFitter evaluation failed with exit code {return_code}")

    likelihood_path = run_dir / "output" / "likelihood.txt"
    if not likelihood_path.is_file():
        raise SystemExit("xFitter did not write output/likelihood.txt")
    likelihood = parse_likelihood(likelihood_path)
    if likelihood["free_parameter_count"] != 0:
        raise SystemExit("Evaluation unexpectedly had free MINUIT parameters")
    if likelihood["nuisance_treatment"] != "fixed":
        raise SystemExit("Evaluation unexpectedly profiled nuisance parameters")
    if likelihood["nuisance_count"] != len(selected_nuisances):
        raise SystemExit("Evaluated nuisance count does not match fixed input")
    projection_pdf = None
    if projection_export is not None:
        projection_path = (
            run_dir
            / "output"
            / str(projection_export["name"])
            / f"{projection_export['name']}_0000.dat"
        )
        if not projection_path.is_file():
            raise SystemExit("xFitter did not write the requested projection PDF grid")
        projection_pdf = {
            "member_file": str(projection_path.relative_to(run_dir)),
            "member_file_sha256": sha256(projection_path),
            "q_ext_GeV": projection_export["q_values"][0],
            "x_point_count": len(projection_export["x_values"]),
            "value_precision": "IEEE-754 double rendered with %.17e",
        }

    group_terms = {}
    for group, dataset_indices in DATASET_GROUPS.items():
        data_chi2 = sum(
            float(likelihood[f"dataset_{index}_data_chi2"])
            for index in dataset_indices
        )
        log_penalty_chi2 = sum(
            float(likelihood[f"dataset_{index}_log_penalty_chi2"])
            for index in dataset_indices
        )
        correlated_penalty_chi2 = sum(
            float(nuisance["shift"]) ** 2
            for nuisance in selected_nuisances
            if nuisance["group"] == group
        )
        group_terms[group] = {
            "data_chi2": data_chi2,
            "correlated_penalty_chi2": correlated_penalty_chi2,
            "log_penalty_chi2": log_penalty_chi2,
            "total_chi2": data_chi2
            + correlated_penalty_chi2
            + log_penalty_chi2,
        }
    grouped_total = sum(group["total_chi2"] for group in group_terms.values())
    if abs(grouped_total - float(likelihood["total_chi2"])) > 1e-9:
        raise SystemExit(
            "HERA+CMS decomposition does not reproduce the combined likelihood: "
            f"{grouped_total} versus {likelihood['total_chi2']}"
        )

    closure = None
    if not overrides:
        expected = float(fit["minimum"]["chi2"])
        scan_config = yaml.safe_load((scan_dir / "scan_config.yaml").read_text())
        tolerance = float(
            scan_config["reference_closure"]["max_absolute_chi2_difference"]
        )
        delta = float(likelihood["total_chi2"]) - expected
        closure = {
            "expected_total_chi2": expected,
            "evaluated_total_chi2": float(likelihood["total_chi2"]),
            "delta_chi2": delta,
            "max_absolute_chi2_difference": tolerance,
            "passed": abs(delta) <= tolerance,
        }

    source_files = [
        script_path,
        project_root / "xfitter" / "include" / "steering.inc",
        project_root / "xfitter" / "src" / "read_steer.f",
        project_root / "xfitter" / "src" / "GetChisquare.f",
        project_root / "xfitter" / "src" / "fcn.f",
        project_root / "xfitter" / "src" / "lhapdf6_output.cc",
    ]
    source_paths = [
        str(path.relative_to(project_root / "xfitter")) for path in source_files[1:]
    ]
    xfitter_patch = subprocess.check_output(
        ["git", "-C", str(project_root / "xfitter"), "diff", "--", *source_paths],
        text=True,
    )
    # The required fixed-nuisance support is committed in the pinned xFitter
    # submodule.  Preserve a local uncommitted diff when present, but a clean
    # worktree is the normal portable-install state and must not be an error.
    xfitter_patch_record = None
    if xfitter_patch.strip():
        patch_path = run_dir / "xfitter_worktree.patch"
        patch_path.write_text(xfitter_patch)
        xfitter_patch_record = {
            "file": patch_path.name,
            "sha256": sha256(patch_path),
        }
    scan_config_path = scan_dir / "scan_config.yaml"
    result = {
        "schema_version": 1,
        "status": "complete" if closure is None or closure["passed"] else "closure_failed",
        "fit_id": fit["fit_id"],
        "parameterization": "direct_HERAPDF",
        "likelihood_group": "combined_with_HERA_CMS_decomposition",
        "minimization": False,
        "profiling": False,
        "parameter_values": values,
        "parameter_overrides": overrides,
        "nuisance_source": str(
            (reference_dir / "nuisances.yaml").relative_to(project_root)
        ),
        "nuisance_stored_precision": nuisance_payload["stored_precision"],
        "nuisance_count": len(selected_nuisances),
        "input_files": input_files,
        "correlation_files": corr_files,
        "likelihood": likelihood,
        "likelihood_groups": group_terms,
        "projection_pdf": projection_pdf,
        "reference_closure": closure,
        "runtime": {
            "started_utc": started.replace(microsecond=0).isoformat(),
            "finished_utc": finished.replace(microsecond=0).isoformat(),
            "wall_seconds": (finished - started).total_seconds(),
        },
        "software": {
            "xfitter_binary": str(xfitter_binary),
            "xfitter_binary_sha256": sha256(xfitter_binary),
            "xfitter_library": str(xfitter_library_path),
            "xfitter_library_sha256": sha256(xfitter_library_path),
            "xfitter_commit": subprocess.check_output(
                ["git", "-C", str(project_root / "xfitter"), "rev-parse", "HEAD"],
                text=True,
            ).strip(),
            "xfitter_worktree_patch": xfitter_patch_record,
            "scan_config_sha256": sha256(scan_config_path),
            "source_sha256": {
                str(path.relative_to(project_root)): sha256(path) for path in source_files
            },
        },
    }
    (run_dir / "evaluation.yaml").write_text(
        yaml.safe_dump(result, sort_keys=False, width=100)
    )
    write_checksums(run_dir)

    print(f"Evaluation stored in {run_dir}")
    print(f"total chi2 = {float(likelihood['total_chi2']):.12f}")
    if closure is not None:
        print(
            "reference closure: "
            f"delta chi2 = {closure['delta_chi2']:+.12g}, "
            f"tolerance = {closure['max_absolute_chi2_difference']:.12g}, "
            f"passed = {closure['passed']}"
        )
        if not closure["passed"]:
            raise SystemExit("Reference closure failed; scan is blocked")


if __name__ == "__main__":
    main()
