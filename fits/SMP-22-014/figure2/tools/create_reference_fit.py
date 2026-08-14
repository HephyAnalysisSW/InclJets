#!/usr/bin/env python3
"""Freeze the completed HERA+CMS global fit into an immutable snapshot."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import yaml


FIT_ID = "smp22014_global_hera_cms_20260805"
FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_revision(repository: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repository), "rev-parse", "HEAD"], text=True
    ).strip()


def parse_parsout(path: Path) -> dict[int, dict[str, float | int | str]]:
    pattern = re.compile(
        r"^\s*(\d+)\s+'([^']+)'\s+"
        r"([+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s+"
        r"([+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s*$"
    )
    result: dict[int, dict[str, float | int | str]] = {}
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            external_id = int(match.group(1))
            result[external_id] = {
                "external_id": external_id,
                "name": match.group(2),
                "value": float(match.group(3)),
                "hesse_error": float(match.group(4)),
            }
    if not result:
        raise RuntimeError(f"No parameters parsed from {path}")
    return result


def parse_covariance(
    path: Path, parameters: dict[int, dict[str, float | int | str]]
) -> tuple[list[dict[str, float | int | str]], np.ndarray]:
    text = path.read_text()
    header = re.search(
        r"EXTERNAL ERROR MATRIX\.\s+NDIM=\s*\d+\s+NPAR=\s*(\d+)", text
    )
    if not header:
        raise RuntimeError("MINUIT covariance header not found")
    size = int(header.group(1))

    correlation_section = text.index("PARAMETER  CORRELATION COEFFICIENTS", header.end())
    matrix_start = text.index("ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.", header.end())
    matrix_start = text.index("\n", matrix_start) + 1
    matrix_text = text[matrix_start:correlation_section]
    values = [float(value) for value in FLOAT_RE.findall(matrix_text)]
    expected = size * (size + 1) // 2
    if len(values) != expected:
        raise RuntimeError(
            f"Expected {expected} lower-triangle values, found {len(values)}"
        )

    correlation_header = re.search(
        r"NO\.[ \t]+GLOBAL[ \t]+((?:\d+[ \t]*)+)",
        text[correlation_section:],
    )
    if not correlation_header:
        raise RuntimeError("MINUIT covariance parameter order not found")
    external_ids = [int(value) for value in correlation_header.group(1).split()]
    if len(external_ids) != size:
        raise RuntimeError(
            f"Expected {size} covariance parameter IDs, found {len(external_ids)}"
        )

    covariance = np.zeros((size, size), dtype=float)
    cursor = 0
    for row in range(size):
        for column in range(row + 1):
            covariance[row, column] = values[cursor]
            covariance[column, row] = values[cursor]
            cursor += 1

    ordered_parameters = [parameters[external_id] for external_id in external_ids]
    return ordered_parameters, covariance


def parse_fit_summary(path: Path, minuit_path: Path) -> dict[str, object]:
    result_text = path.read_text()
    minuit_text = minuit_path.read_text()
    first_line = re.search(
        r"After minimisation\s+([+-]?\d+(?:\.\d+)?)\s+(\d+)\s+([+-]?\d+(?:\.\d+)?)",
        result_text,
    )
    hesse = re.search(
        r"FCN=\s*([+-]?\d+(?:\.\d+)?)\s+FROM HESSE\s+STATUS=(\S+).*?\n"
        r"\s*EDM=\s*([+-]?(?:\d+\.\d*|\.\d+|\d+)(?:E[+-]?\d+))"
        r".*?ERROR MATRIX\s+([^\n]+)",
        minuit_text,
        flags=re.DOTALL,
    )
    if not first_line or not hesse:
        raise RuntimeError("Could not parse fit minimum from Results/MINUIT output")

    datasets = []
    dataset_pattern = re.compile(
        r"^Dataset\s+(\d+)\s+([+-]?\d+\.\d+)\(\s*([+-]?\d+\.\d+)\)"
        r"\s+(\d+)\s+(.+?)\s*$"
    )
    for line in result_text.splitlines():
        match = dataset_pattern.match(line)
        if match:
            datasets.append(
                {
                    "index": int(match.group(1)),
                    "reported_partial_chi2": float(match.group(2)),
                    "reported_shift": float(match.group(3)),
                    "npoints": int(match.group(4)),
                    "name": match.group(5).strip(),
                    "group": "HERA" if int(match.group(1)) <= 7 else "CMS",
                }
            )

    def scalar(pattern: str, cast: type = float) -> object:
        match = re.search(pattern, result_text)
        if not match:
            raise RuntimeError(f"Missing Results.txt field matching {pattern!r}")
        return cast(match.group(1))

    return {
        "chi2": float(hesse.group(1)),
        "ndof": int(first_line.group(2)),
        "chi2_per_dof_reported": float(first_line.group(3)),
        "edm": float(hesse.group(3)),
        "minuit_status": hesse.group(2),
        "covariance_status": hesse.group(4).strip().lower(),
        "errordef": 1.0,
        "reported_correlated_chi2": scalar(r"Correlated Chi2\s+([+-]?\d+(?:\.\d+)?)"),
        "reported_log_penalty_chi2": scalar(
            r"Log penalty Chi2\s+([+-]?\d+(?:\.\d+)?)"
        ),
        "nuisance_count": scalar(r"Systematic shifts\s+(\d+)", int),
        "reported_dataset_terms": datasets,
    }


def parse_nuisances(path: Path) -> list[dict[str, object]]:
    pattern = re.compile(
        r"^\s*(\d+)\s+(.+?)\s+"
        r"([+-]?\d+\.\d+)\s+\+/-\s+([+-]?\d+\.\d+)\s+(\S+)\s*$"
    )
    nuisances = []
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        index = int(match.group(1))
        if not 1 <= index <= 198:
            continue
        nuisances.append(
            {
                "index": index,
                "name": match.group(2).strip(),
                "shift": float(match.group(3)),
                "error": float(match.group(4)),
                "type": match.group(5),
                "group": "HERA" if index <= 162 else "CMS",
            }
        )
    if len(nuisances) != 198:
        raise RuntimeError(f"Expected 198 nuisance parameters, found {len(nuisances)}")
    return nuisances


def copy_files(source: Path, destination: Path, names: list[str]) -> None:
    destination.mkdir(parents=True)
    for name in names:
        shutil.copy2(source / name, destination / name)


def write_yaml(path: Path, payload: object) -> None:
    path.write_text(yaml.safe_dump(payload, sort_keys=False, width=100))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="Snapshot directory (default: figure2/reference_fit)",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    figure2_dir = script_path.parents[1]
    project_root = script_path.parents[4]
    output = (args.output or figure2_dir / "reference_fit").resolve()
    if output.exists():
        raise SystemExit(
            f"Refusing to overwrite immutable snapshot {output}; choose a new --output"
        )

    central_dir = figure2_dir / "central"
    errors_dir = figure2_dir / "errors"
    central_output = central_dir / "output"
    errors_output = errors_dir / "output"

    parameters = parse_parsout(errors_output / "parsout_0")
    ordered_parameters, covariance = parse_covariance(
        errors_output / "minuit.out.txt", parameters
    )
    summary = parse_fit_summary(
        errors_output / "Results.txt", errors_output / "minuit.out.txt"
    )
    nuisances = parse_nuisances(errors_output / "Results.txt")
    if summary["nuisance_count"] != len(nuisances):
        raise RuntimeError("Nuisance count differs between summary and parsed table")

    output.mkdir(parents=True)
    copy_files(
        central_dir,
        output / "cards" / "minimization",
        ["parameters.yaml", "steering.txt", "constants.yaml"],
    )
    copy_files(
        errors_dir,
        output / "cards" / "covariance",
        ["parameters.yaml", "steering.txt", "constants.yaml"],
    )
    copy_files(
        central_output,
        output / "raw" / "minimization",
        ["Results.txt", "minuit.out.txt", "parsout_1", "fittedresults.txt"],
    )
    copy_files(
        errors_output,
        output / "raw" / "covariance",
        ["Results.txt", "minuit.out.txt", "parsout_0", "parsout_1", "fittedresults.txt"],
    )

    parameter_names = np.asarray(
        [str(parameter["name"]) for parameter in ordered_parameters]
    )
    values = np.asarray(
        [float(parameter["value"]) for parameter in ordered_parameters], dtype=float
    )
    errors = np.asarray(
        [float(parameter["hesse_error"]) for parameter in ordered_parameters],
        dtype=float,
    )
    covariance_errors = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(covariance_errors, covariance_errors)
    eigenvalues = np.linalg.eigvalsh(covariance)
    np.savez_compressed(
        output / "covariance.npz",
        parameter_names=parameter_names,
        external_parameter_ids=np.asarray(
            [int(parameter["external_id"]) for parameter in ordered_parameters]
        ),
        values=values,
        hesse_errors=errors,
        covariance=covariance,
        correlation=correlation,
        covariance_eigenvalues=eigenvalues,
    )

    nuisance_payload = {
        "schema_version": 1,
        "fit_id": FIT_ID,
        "description": "Experimental nuisance values at the global HERA+CMS minimum.",
        "source": "raw/covariance/Results.txt",
        "stored_precision": "four decimal places (xFitter Results.txt F9.4 output)",
        "exact_runtime_values_available": False,
        "usage_note": (
            "These are the only retained nuisance values from the completed fit. "
            "Check direct-reference chi2 closure before using them as fixed scan inputs."
        ),
        "count": len(nuisances),
        "group_boundary": {
            "HERA": [1, 162],
            "CMS": [163, 198],
        },
        "nuisances": nuisances,
    }
    write_yaml(output / "nuisances.yaml", nuisance_payload)

    source_files = [
        script_path,
        central_dir / "parameters.yaml",
        central_dir / "steering.txt",
        central_dir / "constants.yaml",
        errors_dir / "parameters.yaml",
        errors_dir / "steering.txt",
        errors_dir / "constants.yaml",
        errors_output / "parsout_0",
        errors_output / "Results.txt",
        errors_output / "minuit.out.txt",
    ]
    source_hashes = {
        str(path.relative_to(project_root)): sha256(path) for path in source_files
    }

    manifest = {
        "schema_version": 1,
        "fit_id": FIT_ID,
        "snapshot_created_utc": dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "status": "immutable_reference",
        "description": "Completed HERA I+II plus CMS inclusive-jet global minimum and Hesse result.",
        "minimum": summary,
        "free_parameters": ordered_parameters,
        "parameter_notes": {
            "count": len(ordered_parameters),
            "independent_pdf_shape_parameters": 15,
            "independent_alpha_s_parameters": 1,
            "dependent_parameters": {
                "Ag": "SUMRULE",
                "Auv": "SUMRULE",
                "Adv": "SUMRULE",
                "AUbar": "=(1-fs)*ADbar",
                "fd": "=1-fs",
            },
            "fixed_parameters": {"fs": 0.4, "ZERO": 0.0},
            "warning": "Dependent values printed as 1.0 in parsout are placeholders; replay the YAML relations.",
        },
        "covariance": {
            "file": "covariance.npz",
            "source": "raw/covariance/minuit.out.txt",
            "parameter_order": parameter_names.tolist(),
            "printed_precision": "MINUIT lower triangle rounded to three significant digits",
            "minimum_eigenvalue_of_stored_rounded_matrix": float(eigenvalues.min()),
            "maximum_eigenvalue_of_stored_rounded_matrix": float(eigenvalues.max()),
            "maximum_relative_diagonal_error_difference": float(
                np.max(np.abs(covariance_errors - errors) / errors)
            ),
        },
        "nuisances": {
            "file": "nuisances.yaml",
            "count": len(nuisances),
            "future_scan_treatment": "fixed_at_global_best_fit; no profiling or minimization",
            "stored_precision": "four decimal places (xFitter Results.txt F9.4 output)",
            "exact_runtime_values_available": False,
            "required_validation": "direct-reference chi2 closure before scanning",
        },
        "theory": {
            "order": "NNLO",
            "evolution": "QCDNUM",
            "q0_GeV": 1.378404875209,
            "n_flavors": 5,
            "flavor_scheme": "variable",
        },
        "software": {
            "xfitter_version": "2.2.1",
            "xfitter_source_commit": git_revision(project_root / "xfitter"),
            "qcdnum_version": "18-00-00",
            "datafiles_commit": git_revision(project_root / "xfitter-datafiles"),
            "platform": "macOS arm64",
        },
        "source_locations": {
            "minimization": "../central",
            "covariance_continuation": "../errors",
            "generator": "../tools/create_reference_fit.py",
        },
        "source_sha256": source_hashes,
    }
    write_yaml(output / "fit_result.yaml", manifest)

    readme = f"""# Immutable global-fit reference

Fit ID: `{FIT_ID}`

This directory freezes the completed HERA I+II plus CMS inclusive-jet global
minimum and its successful 16-parameter Hesse continuation. It is the sole
reference point for the no-minimization likelihood scans.

- `fit_result.yaml`: canonical parameters, fit summary, provenance, and hashes.
- `covariance.npz`: ordered values, Hesse errors, covariance, correlation, and eigenvalues.
- `nuisances.yaml`: all 198 global-best-fit experimental nuisance values as
  retained by xFitter's four-decimal `Results.txt` output.
- `cards/`: exact minimization and covariance-continuation cards.
- `raw/`: minimal xFitter outputs from which the snapshot was extracted.
- `files.sha256`: integrity hashes for every snapshot payload file.

Do not edit or regenerate this directory in place. The generator refuses to
overwrite it; create a new output directory and fit ID for a new reference.
No likelihood evaluation or minimization is performed by the generator.
The exact runtime nuisance array was not saved by the completed fit. Before a
scan uses the rounded values as fixed inputs, its evaluator must reproduce the
direct-reference chi2 within an explicitly recorded tolerance.
"""
    (output / "README.md").write_text(readme)

    payload_files = sorted(
        path for path in output.rglob("*") if path.is_file() and path.name != "files.sha256"
    )
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(output)}" for path in payload_files
    ]
    (output / "files.sha256").write_text("\n".join(checksum_lines) + "\n")

    print(f"Created immutable reference snapshot: {output}")
    print(f"Parameters: {len(ordered_parameters)}; nuisances: {len(nuisances)}")
    print(
        "Stored rounded covariance eigenvalue range: "
        f"[{eigenvalues.min():.6e}, {eigenvalues.max():.6e}]"
    )


if __name__ == "__main__":
    main()
