#!/usr/bin/env python3
"""Scan ADbar with a common Q0=1.65 GeV table route for direct and POD PDFs.

This is a CMS-only, fixed-zero-nuisance *evaluation* scan: no parameter is
minimized or profiled.  It removes the grid/import and Q0-history offsets
before comparing the full-POD representation with the native parametrisation.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml

from run_gluon_isolation import export_full_direct_grid, parameters, write_lhapdf_set


def make_steering(figure_dir: Path, source_direct: Path) -> tuple[str, str]:
    nuisances = yaml.safe_load((figure_dir / "reference_fit" / "nuisances.yaml").read_text())["nuisances"]
    active = [item for item in nuisances if item["group"] == "CMS" and not str(item["name"]).startswith("proc_")]
    if len(active) != 29:
        raise RuntimeError(f"Expected 29 active CMS sources, found {len(active)}")
    fixed = "# local_index source_name fixed_shift\n" + "\n".join(
        f"{index:4d} {item['name']} 0.0" for index, item in enumerate(active, 1)) + "\n"
    steering = (figure_dir / "pod_gluon_cms_7mode" / "steering.txt").read_text()
    steering = steering.replace("&End\n\n&Output", "  UseFixedNuisances = True\n  FixedNuisanceFile = 'fixed_nuisances.dat'\n&End\n\n&Output", 1)
    source = (source_direct / "steering.txt").read_text()
    return steering + "\n" + source[source.index("&Cuts"):], fixed


def run_case(run_dir: Path, source: Path, steering: str, fixed: str, project_root: Path, figure_dir: Path, run_xfitter, parameter_text: str | None = None) -> dict[str, object]:
    run_dir.mkdir(parents=True, exist_ok=True)
    if parameter_text is None:
        shutil.copy2(source / "parameters.yaml", run_dir / "parameters.yaml")
    else:
        (run_dir / "parameters.yaml").write_text(parameter_text)
    for name in ("constants.yaml",):
        shutil.copy2(source / name, run_dir / name)
    (run_dir / "steering.txt").write_text(steering)
    (run_dir / "fixed_nuisances.dat").write_text(fixed)
    (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
    (run_dir / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
    if (source / "pod.yaml").is_file():
        shutil.copy2(source / "pod.yaml", run_dir / "pod.yaml")
    run_xfitter(run_dir, project_root)
    result = parse_likelihood(run_dir / "output" / "likelihood.txt")
    if result["free_parameter_count"] != 0 or result["nuisance_treatment"] != "fixed" or result["nuisance_count"] != 29:
        raise RuntimeError(f"Fixed CMS likelihood contract failed: {result}")
    return result


def points(source: Path) -> list[tuple[str, Path]]:
    candidates = [("reference", source / "runs" / "_reference")]
    candidates += [(path.parent.parent.name, path.parent.parent) for path in (source / "runs" / "ADbar").glob("*/direct/evaluation.yaml")]
    result = []
    for label, root in candidates:
        evaluation = yaml.safe_load((root / "direct" / "evaluation.yaml").read_text())
        result.append((label, root, float(evaluation["parameter_values"]["ADbar"])))
    return sorted(result, key=lambda item: item[2])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Completed production likelihood-scan directory")
    parser.add_argument("--output", type=Path, default=Path("output_scale_aligned_adbar"))
    parser.add_argument("--limit", type=int, help="Run only the first N ordered points (pilot/debugging)")
    args = parser.parse_args()
    script_dir = Path(__file__).resolve().parent
    figure_dir, project_root = script_dir.parent, script_dir.parents[3]
    sys.path.insert(0, str(figure_dir / "likelihood_scans"))
    global parse_likelihood
    from scan_tools import parse_likelihood, run_xfitter
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    todo = points(source)
    if args.limit:
        todo = todo[:args.limit]
    records = []
    old_path = os.environ.get("LHAPDF_DATA_PATH", "")
    try:
        for number, (label, root, adbar) in enumerate(todo, 1):
            point = output / label
            saved = point / "result.yaml"
            if saved.is_file():
                record = yaml.safe_load(saved.read_text()); records.append(record)
                print(f"[{number}/{len(todo)}] {label}: cached")
                continue
            direct, pod = root / "direct", root / "full_pod"
            point.mkdir(parents=True, exist_ok=True)
            steering, fixed = make_steering(figure_dir, direct)
            set_name = f"direct_{label}"
            grid = point / f"direct_export_{set_name}" / "output" / set_name
            if not grid.is_dir():
                grid = export_full_direct_grid(direct, point, project_root, figure_dir, run_xfitter, 1.65, set_name)
            lhapdf_root = point / "lhapdf"
            write_lhapdf_set(grid, lhapdf_root / set_name, set_name)
            os.environ["LHAPDF_DATA_PATH"] = str(lhapdf_root) + (":" + old_path if old_path else "")
            direct_like = run_case(point / "matched_direct", direct, steering, fixed, project_root, figure_dir, run_xfitter, parameters(set_name, 1.65))
            pod_like = run_case(point / "full_pod", pod, steering, fixed, project_root, figure_dir, run_xfitter)
            record = {"label": label, "ADbar": adbar, "matched_direct_chi2": float(direct_like["total_chi2"]), "full_pod_chi2": float(pod_like["total_chi2"]), "pod_minus_direct": float(pod_like["total_chi2"]) - float(direct_like["total_chi2"])}
            saved.write_text(yaml.safe_dump(record, sort_keys=False)); records.append(record)
            print(f"[{number}/{len(todo)}] ADbar={adbar:.7g}: direct={record['matched_direct_chi2']:.6f}, POD={record['full_pod_chi2']:.6f}, delta={record['pod_minus_direct']:+.6f}")
    finally:
        if old_path: os.environ["LHAPDF_DATA_PATH"] = old_path
        else: os.environ.pop("LHAPDF_DATA_PATH", None)
    records.sort(key=lambda item: item["ADbar"])
    (output / "results.yaml").write_text(yaml.safe_dump({"description": __doc__, "source": str(source), "q0_GeV": 1.65, "points": records}, sort_keys=False))
    np.savez(output / "results.npz", ADbar=np.array([x["ADbar"] for x in records]), matched_direct_chi2=np.array([x["matched_direct_chi2"] for x in records]), full_pod_chi2=np.array([x["full_pod_chi2"] for x in records]))


if __name__ == "__main__":
    main()
