#!/usr/bin/env python3
"""Evaluate HERA, CMS, and combined fixed-global-nuisance ADbar terms at Q0=1.65.

The direct input is the already-exported per-point LHAPDF table and the POD
input is the matching full 100-mode reconstruction.  Both use the original
11-dataset steering card and the same stored 198 global-fit nuisance shifts.
No minimization or profiling is performed.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml



def table_parameters(template: Path, set_name: str) -> str:
    """Replace only the PDF input/evolution while retaining DIS reaction setup."""
    text = (template / "parameters.yaml").read_text().split("\nWriteLHAPDF6:", 1)[0].rstrip() + "\n"
    replacement = f"""DefaultDecomposition: proton
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

Order:"""
    text, count = re.subn(r"DefaultDecomposition: proton\nDecompositions:.*?\nOrder:", replacement, text, flags=re.DOTALL)
    if count != 1: raise RuntimeError(f"Could not replace direct PDF setup in {template}")
    return re.sub(r"(?m)^Q0:\s*.*$", "Q0: 1.65", text)


def points(source: Path) -> list[tuple[str, Path, float]]:
    raw = [("reference", source / "runs" / "_reference")]
    raw += [(p.parent.parent.name, p.parent.parent) for p in (source / "runs" / "ADbar").glob("*/direct/evaluation.yaml")]
    return sorted([(label, root, float(yaml.safe_load((root / "direct" / "evaluation.yaml").read_text())["parameter_values"]["ADbar"])) for label, root in raw], key=lambda x: x[2])


def evaluate(target: Path, template: Path, project_root: Path, figure_dir: Path, run_xfitter, parameter_text: str | None = None) -> dict[str, object]:
    target.mkdir(parents=True, exist_ok=True)
    (target / "parameters.yaml").write_text(parameter_text if parameter_text is not None else (template / "parameters.yaml").read_text())
    for name in ("constants.yaml", "steering.txt", "fixed_nuisances.dat"):
        shutil.copy2(template / name, target / name)
    if (template / "pod.yaml").is_file(): shutil.copy2(template / "pod.yaml", target / "pod.yaml")
    (target / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
    (target / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
    run_xfitter(target, project_root)
    return parse_likelihood(target / "output" / "likelihood.txt")


def totals(likelihood: dict[str, object], nuisances: list[dict[str, object]]) -> dict[str, float]:
    grouped = likelihood_groups(likelihood, nuisances)
    return {key: float(value["total_chi2"]) for key, value in grouped.items()} | {"sum": float(likelihood["total_chi2"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--scale-scan", required=True, type=Path, help="Completed scale-aligned CMS scan, providing LHAPDF grids")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent; figure_dir, project_root = here.parent, here.parents[3]
    sys.path.insert(0, str(figure_dir / "likelihood_scans"))
    global parse_likelihood, likelihood_groups
    from scan_tools import parse_likelihood, likelihood_groups, run_xfitter
    source, scale, output = args.source.resolve(), args.scale_scan.resolve(), args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    nuisances = yaml.safe_load((figure_dir / "reference_fit" / "nuisances.yaml").read_text())["nuisances"]
    if len(nuisances) != 198: raise RuntimeError(f"Expected 198 stored nuisances, found {len(nuisances)}")
    old_path, records = os.environ.get("LHAPDF_DATA_PATH", ""), []
    try:
        for index, (label, root, adbar) in enumerate(points(source), 1):
            point, saved = output / label, output / label / "result.yaml"
            if saved.is_file(): records.append(yaml.safe_load(saved.read_text())); print(f"[{index}/21] {label}: cached"); continue
            set_name = f"direct_{label}"
            grid = scale / label / "lhapdf" / set_name
            if not grid.is_dir(): raise RuntimeError(f"Missing table input {grid}")
            os.environ["LHAPDF_DATA_PATH"] = str(grid.parent) + (":" + old_path if old_path else "")
            direct = evaluate(point / "matched_direct", root / "direct", project_root, figure_dir, run_xfitter, table_parameters(root / "direct", set_name))
            pod = evaluate(point / "full_pod", root / "full_pod", project_root, figure_dir, run_xfitter)
            record = {"label": label, "ADbar": adbar, "matched_direct": totals(direct, nuisances), "full_pod": totals(pod, nuisances)}
            saved.write_text(yaml.safe_dump(record, sort_keys=False)); records.append(record)
            print(f"[{index}/21] ADbar={adbar:.7f}: HERA delta={record['full_pod']['HERA']-record['matched_direct']['HERA']:+.3f}, CMS delta={record['full_pod']['CMS']-record['matched_direct']['CMS']:+.3f}")
    finally:
        if old_path: os.environ["LHAPDF_DATA_PATH"] = old_path
        else: os.environ.pop("LHAPDF_DATA_PATH", None)
    records.sort(key=lambda x: x["ADbar"])
    arrays={"ADbar":np.array([x["ADbar"] for x in records])}
    for route in ("matched_direct", "full_pod"):
        for group in ("HERA", "CMS", "sum"): arrays[f"{route}_{group}"]=np.array([x[route][group] for x in records])
    np.savez(output / "results.npz", **arrays)
    (output / "results.yaml").write_text(yaml.safe_dump({"description":__doc__, "points":records}, sort_keys=False))

if __name__ == "__main__": main()
