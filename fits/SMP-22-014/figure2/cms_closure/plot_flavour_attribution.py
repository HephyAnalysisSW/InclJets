#!/usr/bin/env python3
"""Plot flavour-by-flavour CMS POD likelihood attribution from summary.yaml."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import yaml


LABELS = {
    "flavour_+21": r"$g$", "flavour_+2": r"$u$", "flavour_-2": r"$\bar u$",
    "flavour_+1": r"$d$", "flavour_-1": r"$\bar d$", "flavour_+3": r"$s$",
    "flavour_-3": r"$\bar s$", "flavour_+4": r"$c$", "flavour_-4": r"$\bar c$",
    "flavour_+5": r"$b$", "flavour_-5": r"$\bar b$",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    shifts = yaml.safe_load(args.summary.read_text())["total_chi2_shifts"]
    singles = [(LABELS[key], value) for key, value in shifts.items() if key.startswith("flavour_")]
    singles.sort(key=lambda item: item[1])
    labels, values = zip(*singles)
    figure, axis = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
    colors = ["#d62728" if value > 0 else "#1f77b4" for value in values]
    axis.barh(labels, values, color=colors)
    axis.axvline(0, color="0.25", lw=0.8)
    axis.set_xlabel(r"$\Delta\chi^2$: replace one input flavour by full POD")
    axis.set_title("CMS inclusive jets, fixed-zero nuisance shifts")
    sum_single = sum(values)
    all_shift = shifts["all_flavours"]
    axis.text(0.02, 0.02, f"sum singles = {sum_single:+.3f}\nall flavours = {all_shift:+.3f}\ncross term = {all_shift-sum_single:+.3f}", transform=axis.transAxes, va="bottom", fontsize=9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
