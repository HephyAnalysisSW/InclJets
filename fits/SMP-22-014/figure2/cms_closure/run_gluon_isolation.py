#!/usr/bin/env python3
"""Isolate CMS likelihood responses to individual full-POD flavour residuals.

Both evaluations start from an identical LHAPDF input grid at Q=1.65 GeV and
are evolved by QCDNUM.  The hybrid changes only the input gluon: it replaces
selected direct-grid flavours by their full-POD reconstructions.  Hence each
hybrid minus direct-input is a flavour-specific representation effect.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml


CMS_IDS = (4361496, 7992985, 5255963, 1624474)


def parse_rows(path: Path) -> dict[str, np.ndarray]:
    """Extract the four CMS data blocks from xFitter fittedresults.txt."""
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 13:
            continue
        try:
            values = [float(value) for value in fields[:11]]
            dataset_id = int(fields[11])
            plot, pt = fields[12].split("/", 1)
        except ValueError:
            continue
        if dataset_id in CMS_IDS:
            rows.append(values[:11] + [float(dataset_id), float(plot), float(pt)])
    result = np.asarray(rows, dtype=float)
    if result.shape != (78, 14):
        raise RuntimeError(f"Expected 78 CMS rows in {path}; found {result.shape}")
    result = result[np.lexsort((result[:, 13], result[:, 12]))]
    return {
        "pt": result[:, 2], "data": result[:, 3], "uncor": result[:, 4],
        "theory": result[:, 7], "dataset_id": result[:, 11].astype(int),
        "rapidity_bin": result[:, 12].astype(int),
    }


def write_lhapdf_set(source: Path, destination: Path, set_name: str, replacements: dict[int, np.ndarray] | None = None, projected_x: np.ndarray | None = None) -> None:
    """Copy a direct Q=1.65 grid, optionally replacing selected flavours."""
    destination.mkdir(parents=True)
    source_dat = next(source.glob("*.dat"))
    source_info = next(source.glob("*.info"))
    lines = source_dat.read_text().splitlines()
    separator = lines.index("---")
    x_values = np.fromstring(lines[separator + 1], sep=" ")
    q_values = np.fromstring(lines[separator + 2], sep=" ")
    flavours = np.fromstring(lines[separator + 3], sep=" ", dtype=int)
    data_start = separator + 4
    data_stop = data_start + len(q_values) * len(x_values)
    if data_stop > len(lines) or len(flavours) == 0:
        raise RuntimeError("Unexpected LHAPDF grid layout")
    rows = lines[data_start:data_stop]
    if replacements:
        if projected_x is None or any(values.size > x_values.size or not np.allclose(x_values[-values.size:], projected_x) for values in replacements.values()):
            raise RuntimeError("Projected flavour grid is not the high-x part of the exported direct grid")
        columns = {pid: int(np.where(flavours == pid)[0][0]) for pid in replacements}
        for row_index, line in enumerate(rows):
            values = np.fromstring(line, sep=" ")
            if values.shape != flavours.shape:
                raise RuntimeError(f"Malformed LHAPDF row {data_start + row_index + 1}")
            # LHAPDF lhagrid1 stores Q as the innermost index.
            x_index = row_index // len(q_values)
            for pid, replacement in replacements.items():
                if x_index >= len(x_values) - replacement.size:
                    values[columns[pid]] = replacement[x_index - (len(x_values) - replacement.size)]
            rows[row_index] = " ".join(f"{value:.16e}" for value in values)
    if len(q_values) == 1:
        # WriteLHAPDF6 collapses nodes separated by only 1e-6.  The
        # interpolator requires two Q knots even though xFitter asks exactly
        # at Q0, so repeat this input grid at harmless nearby Q values.
        q_values = q_values[0] * np.array([1.0, 1.001, 1.002])
        lines[separator + 2] = " ".join(f"{value:.16e}" for value in q_values)
        rows = [row for row in rows for _ in q_values]
    lines = lines[:data_start] + rows + lines[data_stop:]
    (destination / f"{set_name}_0000.dat").write_text("\n".join(lines) + "\n")
    info = source_info.read_text().replace("Direct HERAPDF scan point", "CMS gluon-isolation direct input")
    if replacements:
        info = info.replace("CMS gluon-isolation direct input", "CMS flavour-selected POD hybrid")
    (destination / f"{set_name}.info").write_text(info)


def parameters(set_name: str, q0: float = 1.65) -> str:
    return f"""Minimizer: MINUIT
MINUIT:
  Commands: |
    call fcn 3
  doErrors: None

Parameters:
  alphas: [0.116638, 0.0]

DefaultDecomposition: proton
Decompositions:
  proton:
    class: LHAPDF
    set: {set_name}
    member: 0

DefaultEvolution: proton-QCDNUM
Evolutions:
  proton-QCDNUM:
    ? !include evolutions/QCDNUM.yaml
    decomposition: proton

Order: NNLO
NFlavour: 5
isFFNS: 0
Q0: {q0:.12g}

? !include constants.yaml

OutputDirectory: output
"""


def export_full_direct_grid(source_direct: Path, output: Path, project_root: Path, figure_dir: Path, run_xfitter, q0: float = 1.65, set_name: str = "cms_direct_full_input") -> Path:
    """Export the direct PDF down to x=1e-9 at a specified input scale."""
    sys.path.insert(0, str(project_root))
    from pod_projection.pod_projection import LHAPDF_XGRID

    run_dir = output / f"direct_export_{set_name}"
    run_dir.mkdir()
    template = (source_direct / "parameters.yaml").read_text()
    template = template.split("\nWriteLHAPDF6:", 1)[0].rstrip()
    template += "\n\n" + yaml.safe_dump({"WriteLHAPDF6": {
        "name": set_name, "evolution": "proton-QCDNUM", "description": "Direct HERAPDF full input grid for CMS controls",
        "Xvalues": [float(value) for value in LHAPDF_XGRID], "Qvalues": [q0, q0 * 1.001, q0 * 1.002],
    }}, sort_keys=False, width=120)
    (run_dir / "parameters.yaml").write_text(template)
    for filename in ("constants.yaml", "steering.txt", "fixed_nuisances.dat"):
        shutil.copy2(source_direct / filename, run_dir / filename)
    (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
    (run_dir / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
    run_xfitter(run_dir, project_root)
    grid = run_dir / "output" / set_name
    if not grid.is_dir():
        raise RuntimeError(f"Missing full direct LHAPDF export: {grid}")
    return grid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Completed full-POD likelihood-scan output")
    parser.add_argument("--output", type=Path, default=Path("output_gluon_isolation"))
    parser.add_argument("--all-flavors", action="store_true", help="Evaluate every single flavour and the all-flavour hybrid")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    figure_dir, project_root = script_dir.parent, script_dir.parents[3]
    scan_dir = figure_dir / "likelihood_scans"
    sys.path.insert(0, str(scan_dir))
    from scan_tools import parse_likelihood, run_xfitter

    source = args.source.resolve()
    source_direct = source / "runs" / "_reference" / "direct"
    projection_path = source / "runs" / "_reference" / "full_pod_projection.npz"
    if not projection_path.is_file():
        raise SystemExit("Source needs full_pod_projection.npz")
    with np.load(projection_path, allow_pickle=False) as projection:
        flavours, x, projected = projection["flavors"].astype(int), projection["x_grid"], projection["projected_grid"]
    projected_by_flavour = {int(pid): projected[index] for index, pid in enumerate(flavours)}

    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)
    direct_grid = export_full_direct_grid(source_direct, output, project_root, figure_dir, run_xfitter)
    lhapdf_root = output / "lhapdf"
    cases: list[tuple[str, dict[int, np.ndarray]]] = [("direct_input", {})]
    if args.all_flavors:
        cases.extend((f"flavour_{pid:+d}", {int(pid): projected_by_flavour[int(pid)]}) for pid in flavours)
        cases.append(("all_flavours", projected_by_flavour))
    else:
        cases.append(("gluon_hybrid", {21: projected_by_flavour[21]}))
    for label, replacements in cases:
        set_name = f"cms_{label}"
        write_lhapdf_set(direct_grid, lhapdf_root / set_name, set_name, replacements or None, x)

    combined = yaml.safe_load((figure_dir / "reference_fit" / "nuisances.yaml").read_text())["nuisances"]
    active_nuisances = [item for item in combined if item["group"] == "CMS" and not str(item["name"]).startswith("proc_")]
    if len(active_nuisances) != 29:
        raise RuntimeError(f"Expected 29 active CMS-only nuisances, found {len(active_nuisances)}")
    fixed = "# local_index source_name fixed_shift\n" + "\n".join(
        f"{index:4d} {item['name']} 0.0" for index, item in enumerate(active_nuisances, start=1)
    ) + "\n"
    steering = (figure_dir / "pod_gluon_cms_7mode" / "steering.txt").read_text()
    steering = steering.replace("&End\n\n&Output", "  UseFixedNuisances = True\n  FixedNuisanceFile = 'fixed_nuisances.dat'\n&End\n\n&Output", 1)
    full_steering = (source_direct / "steering.txt").read_text()
    steering += "\n" + full_steering[full_steering.index("&Cuts"):]

    old_lhapdf_path = os.environ.get("LHAPDF_DATA_PATH", "")
    os.environ["LHAPDF_DATA_PATH"] = str(lhapdf_root) + (":" + old_lhapdf_path if old_lhapdf_path else "")
    records, likelihoods = {}, {}
    try:
        for label, _ in cases:
            set_name = f"cms_{label}"
            run_dir = output / label
            run_dir.mkdir()
            (run_dir / "parameters.yaml").write_text(parameters(set_name))
            shutil.copy2(source_direct / "constants.yaml", run_dir / "constants.yaml")
            (run_dir / "steering.txt").write_text(steering)
            (run_dir / "fixed_nuisances.dat").write_text(fixed)
            (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
            (run_dir / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
            run_xfitter(run_dir, project_root)
            likelihood = parse_likelihood(run_dir / "output" / "likelihood.txt")
            if likelihood["free_parameter_count"] != 0 or likelihood["nuisance_treatment"] != "fixed" or likelihood["nuisance_count"] != 29:
                raise RuntimeError(f"Fixed CMS contract failed for {label}: {likelihood}")
            records[label], likelihoods[label] = parse_rows(run_dir / "output" / "fittedresults.txt"), likelihood
    finally:
        if old_lhapdf_path:
            os.environ["LHAPDF_DATA_PATH"] = old_lhapdf_path
        else:
            os.environ.pop("LHAPDF_DATA_PATH", None)

    direct = records["direct_input"]
    for label, record in records.items():
        if not np.array_equal(direct["pt"], record["pt"]):
            raise RuntimeError(f"Direct-input and {label} CMS rows do not align")
    np.savez_compressed(output / "flavour_isolation.npz", **{f"{label}_{key}": value for label, record in records.items() for key, value in record.items()})
    shifts = {label: float(likelihood["total_chi2"]) - float(likelihoods["direct_input"]["total_chi2"]) for label, likelihood in likelihoods.items() if label != "direct_input"}
    (output / "summary.yaml").write_text(yaml.safe_dump({
        "source": str(source), "method": "same Q=1.65-GeV full direct input grid and QCDNUM evolution; replace selected flavour(s) by full POD",
        "nuisance_treatment": "29 active CMS-only sources fixed to zero", "direct_input_likelihood": likelihoods["direct_input"],
        "hybrid_likelihoods": {label: likelihood for label, likelihood in likelihoods.items() if label != "direct_input"},
        "total_chi2_shifts": shifts,
    }, sort_keys=False, width=110))
    for label, shift in shifts.items():
        print(f"{label} CMS chi2 shift: {shift:.12g}")
    print(f"Wrote {output / 'flavour_isolation.npz'}")


if __name__ == "__main__":
    main()
