#!/usr/bin/env python3
"""Separate native-PDF, grid-import, and input-scale effects in CMS closure."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import yaml

from run_gluon_isolation import export_full_direct_grid, parameters, write_lhapdf_set


def native_q0(parameters_path: Path) -> float:
    match = re.search(r"(?m)^Q0:\s*([0-9.eE+-]+)\s*$", parameters_path.read_text())
    if not match:
        raise RuntimeError(f"Could not find Q0 in {parameters_path}")
    return float(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Completed full-POD likelihood-scan output")
    parser.add_argument("--output", type=Path, default=Path("output_technical_controls"))
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    figure_dir, project_root = script_dir.parent, script_dir.parents[3]
    scan_dir = figure_dir / "likelihood_scans"
    sys.path.insert(0, str(scan_dir))
    from scan_tools import parse_likelihood, run_xfitter

    source = args.source.resolve()
    source_direct = source / "runs" / "_reference" / "direct"
    if not (source_direct / "parameters.yaml").is_file():
        raise SystemExit(f"Missing direct reference: {source_direct}")
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)
    q_native, q_rebased = native_q0(source_direct / "parameters.yaml"), 1.65

    combined = yaml.safe_load((figure_dir / "reference_fit" / "nuisances.yaml").read_text())["nuisances"]
    active = [item for item in combined if item["group"] == "CMS" and not str(item["name"]).startswith("proc_")]
    if len(active) != 29:
        raise RuntimeError(f"Expected 29 active CMS sources, found {len(active)}")
    fixed = "# local_index source_name fixed_shift\n" + "\n".join(f"{i:4d} {item['name']} 0.0" for i, item in enumerate(active, 1)) + "\n"
    steering = (figure_dir / "pod_gluon_cms_7mode" / "steering.txt").read_text()
    steering = steering.replace("&End\n\n&Output", "  UseFixedNuisances = True\n  FixedNuisanceFile = 'fixed_nuisances.dat'\n&End\n\n&Output", 1)
    source_steering = (source_direct / "steering.txt").read_text()
    steering += "\n" + source_steering[source_steering.index("&Cuts"):]

    # These exports use the native analytic PDF; the following table runs are
    # the only locations where LHAPDF interpolation is introduced.
    grid_native = export_full_direct_grid(source_direct, output, project_root, figure_dir, run_xfitter, q_native, "cms_direct_qnative")
    grid_rebased = export_full_direct_grid(source_direct, output, project_root, figure_dir, run_xfitter, q_rebased, "cms_direct_q165")
    lhapdf_root = output / "lhapdf"
    write_lhapdf_set(grid_native, lhapdf_root / "cms_direct_qnative", "cms_direct_qnative")
    write_lhapdf_set(grid_rebased, lhapdf_root / "cms_direct_q165", "cms_direct_q165")

    cases = [("native", None, q_native), ("table_native_q0", "cms_direct_qnative", q_native), ("table_q165", "cms_direct_q165", q_rebased)]
    old_lhapdf_path = os.environ.get("LHAPDF_DATA_PATH", "")
    os.environ["LHAPDF_DATA_PATH"] = str(lhapdf_root) + (":" + old_lhapdf_path if old_lhapdf_path else "")
    likelihoods = {}
    try:
        for label, set_name, q0 in cases:
            run_dir = output / label
            run_dir.mkdir()
            if set_name is None:
                shutil.copy2(source_direct / "parameters.yaml", run_dir / "parameters.yaml")
            else:
                (run_dir / "parameters.yaml").write_text(parameters(set_name, q0))
            shutil.copy2(source_direct / "constants.yaml", run_dir / "constants.yaml")
            (run_dir / "steering.txt").write_text(steering)
            (run_dir / "fixed_nuisances.dat").write_text(fixed)
            (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
            (run_dir / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
            run_xfitter(run_dir, project_root)
            likelihood = parse_likelihood(run_dir / "output" / "likelihood.txt")
            if likelihood["free_parameter_count"] != 0 or likelihood["nuisance_treatment"] != "fixed" or likelihood["nuisance_count"] != 29:
                raise RuntimeError(f"Fixed CMS contract failed for {label}: {likelihood}")
            likelihoods[label] = likelihood
    finally:
        if old_lhapdf_path:
            os.environ["LHAPDF_DATA_PATH"] = old_lhapdf_path
        else:
            os.environ.pop("LHAPDF_DATA_PATH", None)

    chi2 = {label: float(value["total_chi2"]) for label, value in likelihoods.items()}
    summary = {
        "source": str(source), "native_q0_GeV": q_native, "rebased_q0_GeV": q_rebased,
        "nuisance_treatment": "29 active CMS-only sources fixed to zero", "likelihoods": likelihoods,
        "native_to_table_at_native_q0": chi2["table_native_q0"] - chi2["native"],
        "table_native_q0_to_table_q165": chi2["table_q165"] - chi2["table_native_q0"],
        "native_to_table_q165": chi2["table_q165"] - chi2["native"],
    }
    (output / "summary.yaml").write_text(yaml.safe_dump(summary, sort_keys=False, width=110))
    print(yaml.safe_dump({key: summary[key] for key in ("native_to_table_at_native_q0", "table_native_q0_to_table_q165", "native_to_table_q165")}, sort_keys=False))


if __name__ == "__main__":
    main()
