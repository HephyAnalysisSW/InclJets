#!/usr/bin/env python3
"""Evaluate direct and full-POD CMS jets with every nuisance fixed to zero.

The input reference directory must be a completed full-POD likelihood-scan
reference (``runs/_reference``).  Its direct parameters and full-POD
coefficients are reused verbatim, so this is an interface/representation
closure test, not a fit or a projection step.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml


CMS_IDS = (4361496, 7992985, 5255963, 1624474)
CMS_LABELS = (r"$|y|<0.5$", r"$0.5<|y|<1.0$", r"$1.0<|y|<1.5$", r"$1.5<|y|<2.0$")


def parse_rows(path: Path) -> dict[str, np.ndarray]:
    """Extract the four CMS data blocks from xFitter's fittedresults.txt."""
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
        if dataset_id not in CMS_IDS:
            continue
        # ylow, yhigh, pTlow, data, uncor, total, theory original/modified,
        # theory errors, pull, rapidity-plot index, pT from the plot key.
        rows.append(values[:11] + [float(dataset_id), float(plot), float(pt)])
    result = np.asarray(rows, dtype=float)
    if result.shape != (78, 14):
        raise RuntimeError(f"Expected 78 CMS rows in {path}; found {result.shape}")
    result = result[np.lexsort((result[:, 13], result[:, 12]))]
    return {
        "ylow": result[:, 0], "yhigh": result[:, 1], "pt": result[:, 2],
        "data": result[:, 3], "uncor": result[:, 4], "total": result[:, 5],
        "theory_original": result[:, 6], "theory": result[:, 7], "pull": result[:, 10],
        "dataset_id": result[:, 11].astype(int), "rapidity_bin": result[:, 12].astype(int),
    }


def write_zero_nuisances(path: Path, nuisances: list[dict[str, object]]) -> None:
    lines = ["# local_index source_name fixed_shift", "# all CMS correlated shifts fixed to zero"]
    lines.extend(f"{index:4d} {item['name']} 0.0" for index, item in enumerate(nuisances, start=1))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Completed likelihood-scan output directory")
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    figure_dir = script_dir.parent
    project_root = script_dir.parents[3]
    scan_dir = figure_dir / "likelihood_scans"
    sys.path.insert(0, str(scan_dir))
    from scan_tools import parse_likelihood, run_xfitter

    source = args.source.resolve()
    source_direct, source_pod = source / "runs" / "_reference" / "direct", source / "runs" / "_reference" / "full_pod"
    for required in (source_direct / "parameters.yaml", source_pod / "parameters.yaml", source_pod / "pod.yaml"):
        if not required.is_file():
            raise SystemExit(f"Missing completed scan reference input: {required}")
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)

    combined_cms_nuisances = [item for item in yaml.safe_load((figure_dir / "reference_fit" / "nuisances.yaml").read_text())["nuisances"] if item["group"] == "CMS"]
    # The combined reference records 36 CMS-labelled shifts.  Its seven
    # ``proc_*`` entries are combined-fit process terms, not active sources in
    # the CMS-only card.  The remaining 29 data-systematic sources are fixed.
    cms_nuisances = [item for item in combined_cms_nuisances if not str(item["name"]).startswith("proc_")]
    if len(cms_nuisances) != 29:
        raise RuntimeError(f"Expected 29 active CMS-only nuisances, found {len(cms_nuisances)}")
    steering = (figure_dir / "pod_gluon_cms_7mode" / "steering.txt").read_text()
    steering = steering.replace("&End\n\n&Output", "  UseFixedNuisances = True\n  FixedNuisanceFile = 'fixed_nuisances.dat'\n&End\n\n&Output", 1)
    # xFitter's jet reader still opens the Cuts namelist, even though none of
    # these HERA process cuts selects a CMS jet point.  Reuse the benign
    # Figure-2 cuts block so the namelist is complete.
    full_steering = (source_direct / "steering.txt").read_text()
    steering += "\n" + full_steering[full_steering.index("&Cuts"):]

    records: dict[str, dict[str, np.ndarray]] = {}
    likelihoods: dict[str, dict[str, object]] = {}
    for label, source_dir in (("direct", source_direct), ("pod", source_pod)):
        run_dir = output / label
        run_dir.mkdir()
        shutil.copy2(source_dir / "parameters.yaml", run_dir / "parameters.yaml")
        if label == "pod":
            shutil.copy2(source_dir / "pod.yaml", run_dir / "pod.yaml")
        shutil.copy2(source_dir / "constants.yaml", run_dir / "constants.yaml")
        (run_dir / "steering.txt").write_text(steering)
        write_zero_nuisances(run_dir / "fixed_nuisances.dat", cms_nuisances)
        (run_dir / "datafiles").symlink_to(project_root / "xfitter-datafiles", target_is_directory=True)
        (run_dir / "unpolarised.wgt").symlink_to(figure_dir / "central" / "unpolarised.wgt")
        run_xfitter(run_dir, project_root)
        likelihood = parse_likelihood(run_dir / "output" / "likelihood.txt")
        if likelihood["free_parameter_count"] != 0 or likelihood["nuisance_treatment"] != "fixed" or likelihood["nuisance_count"] != 29:
            raise RuntimeError(f"CMS zero-systematics contract failed for {label}: {likelihood}")
        records[label] = parse_rows(run_dir / "output" / "fittedresults.txt")
        likelihoods[label] = likelihood

    if not np.array_equal(records["direct"]["dataset_id"], records["pod"]["dataset_id"]) or not np.array_equal(records["direct"]["pt"], records["pod"]["pt"]):
        raise RuntimeError("Direct and POD CMS rows do not align")
    direct, pod = records["direct"], records["pod"]
    np.savez_compressed(
        output / "cms_zero_systematics.npz",
        **{f"direct_{key}": value for key, value in direct.items()},
        **{f"pod_{key}": value for key, value in pod.items()},
        direct_total_chi2=float(likelihoods["direct"]["total_chi2"]),
        pod_total_chi2=float(likelihoods["pod"]["total_chi2"]),
        direct_data_chi2=float(likelihoods["direct"]["data_chi2"]),
        pod_data_chi2=float(likelihoods["pod"]["data_chi2"]),
        direct_log_penalty_chi2=float(likelihoods["direct"]["log_penalty_chi2"]),
        pod_log_penalty_chi2=float(likelihoods["pod"]["log_penalty_chi2"]),
    )
    (output / "metadata.yaml").write_text(yaml.safe_dump({
        "source": str(source), "datasets": list(CMS_IDS), "nuisance_treatment": "all 29 CMS-only active nuisance shifts fixed to zero",
        "minimization": False, "profiling": False, "direct_likelihood": likelihoods["direct"], "pod_likelihood": likelihoods["pod"],
    }, sort_keys=False, width=110))
    print(f"Wrote {output / 'cms_zero_systematics.npz'}")


if __name__ == "__main__":
    main()
