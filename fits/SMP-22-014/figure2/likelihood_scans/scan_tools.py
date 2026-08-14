#!/usr/bin/env python3
"""Shared, side-effect-free helpers for fixed-nuisance likelihood scans."""

from __future__ import annotations

import hashlib
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import numpy as np


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


def xfitter_environment(project_root: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = ":".join(
        [
            str(project_root / "install" / "xfitter" / "bin"),
            str(project_root / "install" / "qcdnum" / "bin"),
            env.get("PATH", ""),
        ]
    )
    library_variable = "DYLD_LIBRARY_PATH" if platform.system() == "Darwin" else "LD_LIBRARY_PATH"
    external_lib = f"{os.environ['CONDA_PREFIX']}/lib" if "CONDA_PREFIX" in os.environ else ""
    env[library_variable] = ":".join(
        [
            str(project_root / "install" / "xfitter" / "lib"),
            str(project_root / "install" / "qcdnum" / "lib"),
            external_lib,
            env.get(library_variable, ""),
        ]
    )
    return env


def xfitter_library(project_root: Path) -> Path:
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    return project_root / "install" / "xfitter" / "lib" / f"libxfitter{suffix}"


def run_xfitter(run_dir: Path, project_root: Path) -> None:
    binary = project_root / "install" / "xfitter" / "bin" / "xfitter"
    with (run_dir / "run.log").open("w") as log:
        process = subprocess.Popen(
            [str(binary)],
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
    if return_code != 0:
        raise RuntimeError(f"xFitter failed with exit code {return_code}")


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
        elif key in float_fields or (
            key.startswith("dataset_")
            and key.endswith(("_data_chi2", "_log_penalty_chi2"))
        ):
            result[key] = float(value)
        else:
            result[key] = value
    missing = (integer_fields | float_fields | {"nuisance_treatment"}) - result.keys()
    if missing:
        raise RuntimeError(f"Missing likelihood fields: {sorted(missing)}")
    return result


def likelihood_groups(
    likelihood: dict[str, object], nuisances: list[dict[str, object]]
) -> dict[str, dict[str, float]]:
    terms: dict[str, dict[str, float]] = {}
    for group, indices in DATASET_GROUPS.items():
        data = sum(float(likelihood[f"dataset_{i}_data_chi2"]) for i in indices)
        log = sum(
            float(likelihood[f"dataset_{i}_log_penalty_chi2"]) for i in indices
        )
        correlated = sum(
            float(nuisance["shift"]) ** 2
            for nuisance in nuisances
            if nuisance["group"] == group
        )
        terms[group] = {
            "data_chi2": data,
            "correlated_penalty_chi2": correlated,
            "log_penalty_chi2": log,
            "total_chi2": data + correlated + log,
        }
    grouped = sum(term["total_chi2"] for term in terms.values())
    if abs(grouped - float(likelihood["total_chi2"])) > 1.0e-9:
        raise RuntimeError(
            f"HERA+CMS terms ({grouped}) do not close to total ({likelihood['total_chi2']})"
        )
    return terms


def write_fixed_nuisances(
    path: Path,
    nuisances: list[dict[str, object]],
    fit_id: str,
    stored_precision: str,
) -> None:
    lines = [
        "# local_index source_name fixed_shift",
        f"# source_fit_id {fit_id}",
        f"# stored_precision {stored_precision}",
    ]
    lines.extend(
        f"{index:4d} {nuisance['name']} {float(nuisance['shift']):.17g}"
        for index, nuisance in enumerate(nuisances, start=1)
    )
    path.write_text("\n".join(lines) + "\n")


def read_lhagrid_first_q(
    path: Path,
    requested_flavors: Sequence[int],
    expected_x: Sequence[float],
    expected_q: float,
) -> np.ndarray:
    """Read exact x*f values at the first Q node of a one-subgrid LHAPDF file."""
    lines = path.read_text().splitlines()
    try:
        marker = lines.index("---")
    except ValueError as error:
        raise RuntimeError(f"No LHAPDF grid marker in {path}") from error
    x_grid = np.fromstring(lines[marker + 1], sep=" ")
    q_grid = np.fromstring(lines[marker + 2], sep=" ")
    flavors = np.fromstring(lines[marker + 3], sep=" ", dtype=int)
    end = lines.index("---", marker + 1)
    rows = np.asarray(
        [np.fromstring(line, sep=" ") for line in lines[marker + 4 : end]],
        dtype=float,
    )
    expected_x = np.asarray(expected_x, dtype=float)
    if x_grid.shape != expected_x.shape or not np.array_equal(x_grid, expected_x):
        raise RuntimeError("Exported LHAPDF x nodes do not exactly match the POD grid")
    if q_grid.size < 1 or q_grid[0] != expected_q:
        raise RuntimeError(
            f"First exported Q node is {q_grid[0] if q_grid.size else None}, expected {expected_q}"
        )
    expected_rows = x_grid.size * q_grid.size
    if rows.shape != (expected_rows, flavors.size):
        raise RuntimeError(
            f"Unexpected LHAPDF payload shape {rows.shape}; expected {(expected_rows, flavors.size)}"
        )
    cube = rows.reshape(x_grid.size, q_grid.size, flavors.size)
    flavor_index = {int(pid): i for i, pid in enumerate(flavors)}
    missing = [pid for pid in requested_flavors if pid not in flavor_index]
    if missing:
        raise RuntimeError(f"Exported grid lacks requested flavors {missing}")
    return np.asarray(
        [cube[:, 0, flavor_index[int(pid)]] for pid in requested_flavors],
        dtype=float,
    )


def write_checksums(run_dir: Path) -> None:
    payloads = sorted(
        path
        for path in run_dir.rglob("*")
        if path.is_file() and not path.is_symlink() and path.name != "files.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(run_dir)}" for path in payloads]
    (run_dir / "files.sha256").write_text("\n".join(lines) + "\n")
