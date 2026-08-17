#!/usr/bin/env python3
"""Outer 16-parameter analytic fit evaluated through a full POD projection."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import minimize


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("nominal_full_pod"))
    p.add_argument("--maxiter", type=int, default=200)
    p.add_argument("--step", type=float, default=0.05, help="Initial simplex step in Hesse sigma.")
    p.add_argument("--prepare-only", action="store_true")
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    fig, scans, root = here.parent, here.parent / "likelihood_scans", here.parents[3]
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    evals = out / "evaluations"; evals.mkdir(exist_ok=True)
    history_file = out / "history.json"
    history = json.loads(history_file.read_text()) if history_file.exists() else []
    sys.path[:0] = [str(root), str(scans)]
    from pod_projection.pod_projection import LHAPDF_XGRID
    from projection_metrics import Figure2ProjectionOperator
    from scan_tools import read_lhagrid_first_q

    ref = yaml.safe_load((fig / "reference_fit" / "fit_result.yaml").read_text())
    cfg = yaml.safe_load((scans / "scan_config.yaml").read_text())["projection"]
    names = [x["name"] for x in ref["free_parameters"]]
    vals = np.array([float(x["value"]) for x in ref["free_parameters"]])
    errs = np.array([float(x["hesse_error"]) for x in ref["free_parameters"]])
    lo, hi = cfg["x_slice"]; x = LHAPDF_XGRID[lo:hi]; flavors = tuple(map(int, cfg["flavors"]))
    op = Figure2ProjectionOperator.build(
        cfg["basis_set"], int(cfg["coefficient_count"]), flavors, x, float(cfg["q_ext_GeV"]),
        cfg["metric"], relative_weight=cfg["relative_weight"],
        relative_x_range=tuple(cfg["relative_x_range"]), relative_floor=cfg["relative_floor"],
        relative_valence_weight=cfg["relative_valence_weight"],
        relative_valence_x_range=tuple(cfg["relative_valence_x_range"]),
        relative_valence_floor=cfg["relative_valence_floor"],
        relative_f2_weight=cfg["relative_f2_weight"],
        relative_f2_x_range=tuple(cfg["relative_f2_x_range"]), relative_f2_floor=cfg["relative_f2_floor"],
    )

    def objective(z: np.ndarray) -> float:
        physical = vals + np.asarray(z) * errs
        key = ",".join(f"{v:.12g}" for v in physical)
        for row in history:
            if row["key"] == key: return float(row["chi2"])
        index = len(history) + 1; point = evals / f"eval_{index:05d}"; direct = point / "direct"; pod = point / "pod"; point.mkdir()
        direct_cmd = [sys.executable, str(scans / "evaluate_direct.py"), "--export-projection-pdf", "--output", str(direct)]
        for name, value in zip(names, physical): direct_cmd += ["--parameter", f"{name}={value:.17g}"]
        with (point / "direct.console.log").open("w") as log:
            subprocess.run(direct_cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        member = direct / "output" / "direct_projection_target" / "direct_projection_target_0000.dat"
        target = read_lhagrid_first_q(member, flavors, x, float(cfg["q_ext_GeV"]))
        projected, coeff, residual = op.project_grid(target)
        coeff_file = point / "coefficients.npz"
        np.savez_compressed(coeff_file, coefficients=coeff, target_grid=target, projected_grid=projected, residual_grid=residual, flavors=flavors, x_grid=x)
        pod_cmd = [sys.executable, str(scans / "evaluate_pod.py"), "--profile-nuisances", "--coefficients", str(coeff_file), "--alphas", f"{physical[names.index('alphas')]:.17g}", "--source-direct-evaluation", str(direct), "--output", str(pod)]
        with (point / "pod.console.log").open("w") as log:
            subprocess.run(pod_cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        result = yaml.safe_load((pod / "evaluation.yaml").read_text()); chi2 = float(result["likelihood"]["total_chi2"])
        history.append({"index": index, "key": key, "z": list(map(float, z)), "parameters": dict(zip(names, map(float, physical))), "chi2": chi2, "directory": str(point.relative_to(out))})
        history_file.write_text(json.dumps(history, indent=2) + "\n")
        print(f"[{index}] chi2={chi2:.8f}; best={min(float(r['chi2']) for r in history):.8f}", flush=True)
        return chi2

    z0 = np.asarray(min(history, key=lambda r: float(r["chi2"]))["z"]) if history else np.zeros(len(names))
    initial = objective(z0)
    if a.prepare_only:
        print(f"Prepared central point, chi2={initial:.8f}"); return
    simplex = np.vstack([z0, *(z0 + np.eye(len(names))[i] * a.step for i in range(len(names)))])
    started = time.time()
    result = minimize(objective, z0, method="Nelder-Mead", options={"maxiter": a.maxiter, "xatol": 1e-3, "fatol": 1e-3, "initial_simplex": simplex, "disp": True})
    summary = {
        "status": "complete" if result.success else "stopped", "message": str(result.message),
        "outer_method": "Nelder-Mead", "evaluation_count": len(history), "wall_seconds": time.time()-started,
        "reference_fit_id": ref["fit_id"], "reference_chi2": float(ref["minimum"]["chi2"]),
        "pod_mediated_chi2": float(result.fun), "z_best": list(map(float, result.x)),
        "parameters_best": dict(zip(names, map(float, vals + result.x * errs))),
        "profiling": "xFitter Hessian nuisance profiling at each POD likelihood point", "projection": cfg,
    }
    (out / "fit_summary.yaml").write_text(yaml.safe_dump(summary, sort_keys=False))
    print(yaml.safe_dump(summary, sort_keys=False))


if __name__ == "__main__":
    main()
