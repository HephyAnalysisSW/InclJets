#!/usr/bin/env python3
"""Plot the scale-aligned CMS ADbar scan in the standard likelihood-scan style."""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
with np.load(args.input, allow_pickle=False) as data:
    x, direct, pod = data["ADbar"], data["matched_direct_chi2"], data["full_pod_chi2"]
zero = int(np.argmin(np.abs(x - 0.267003)))
sigma = (x - x[zero]) / ((x[-1] - x[0]) / 4.0)
fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True, constrained_layout=True, gridspec_kw={"height_ratios": [2.2, 1]})
baseline = direct[zero]
axes[0].plot(sigma, direct-baseline, "o-", label="matched direct table, $Q_0=1.65$ GeV", color="#1f77b4")
axes[0].plot(sigma, pod-baseline, "s--", label="full POD, $Q_0=1.65$ GeV", color="#d62728")
axes[0].axhline(0, color=".72", lw=.8); axes[0].axvline(0, color=".82", lw=.8)
axes[0].set_ylabel(r"$\chi^2-\chi^2_{\mathrm{direct}}(0)$"); axes[0].grid(alpha=.25); axes[0].legend(frameon=False, fontsize=9)
axes[1].plot(sigma, pod-direct, "D-", color="#6a3d9a")
axes[1].axhline(0, color=".35", lw=.8); axes[1].axvline(0, color=".82", lw=.8)
axes[1].set_xlabel(r"ADbar scan coordinate [$\sigma_{\mathrm{HESSE}}$]"); axes[1].set_ylabel(r"$\chi^2_{\rm POD}-\chi^2_{\rm direct}$"); axes[1].grid(alpha=.25)
fig.suptitle("CMS inclusive jets: scale-aligned full-POD closure scan\nfixed-zero-nuisance likelihood; reference ADbar = 0.267003")
args.output.parent.mkdir(parents=True, exist_ok=True); fig.savefig(args.output, dpi=180)
print(f"Wrote {args.output}")
