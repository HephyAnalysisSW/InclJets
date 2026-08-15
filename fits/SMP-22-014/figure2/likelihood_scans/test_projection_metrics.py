#!/usr/bin/env python3
"""Deterministic unit check for the Figure-2 POD projection objectives.

This test deliberately uses a synthetic POD matrix.  It verifies the
Figure-2-only relative gluon, valence and F2 equations without LHAPDF grids,
xFitter, or public data downloads.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from projection_metrics import Figure2ProjectionOperator


def main() -> None:
    flavors = (21, 2, -2, 1, -1, 3, -3, 4, -4, 5, -5)
    x_grid = np.array((1.0e-4, 5.0e-3, 0.06, 0.4, 0.95))
    n_modes = 3
    rng = np.random.default_rng(22014)
    matrix = rng.normal(size=(len(flavors) * len(x_grid), n_modes))
    reference = rng.normal(size=(len(flavors), len(x_grid)))
    base = SimpleNamespace(
        basis=SimpleNamespace(flavors=flavors),
        x_grid=x_grid,
        matrix=matrix,
        reference_grid=reference,
    )
    operator = Figure2ProjectionOperator(
        base=base,
        metric="relative_gluon",
        relative_weight=0.1,
        relative_x_range=(0.05, 0.99),
        relative_valence_weight=0.25,
        relative_valence_x_range=(1.0e-4, 0.1),
        relative_f2_weight=3.0,
        relative_f2_x_range=(1.0e-4, 0.1),
    )
    expected = np.array((0.17, -0.32, 0.41))
    target = reference + (matrix @ expected).reshape(reference.shape)
    projected, coefficients, residual = operator.project_grid(target)

    np.testing.assert_allclose(coefficients, expected, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(projected, target, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(residual, 0.0, rtol=0.0, atol=2.0e-14)
    print("Figure-2 relative POD projection metric: PASS")


if __name__ == "__main__":
    main()
