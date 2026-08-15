"""Figure-2-specific closure objectives for POD likelihood scans.

These terms change only the projection of an external PDF onto the POD basis.
They are neither experimental likelihood terms nor PDF priors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from pod_projection.pod_projection import ProjectionOperator


@dataclass(frozen=True)
class Figure2ProjectionOperator:
    """POD operator augmented with the Figure-2 PDF-closure objective."""

    base: ProjectionOperator
    metric: str
    relative_weight: float = 0.1
    relative_x_range: tuple[float, float] = (0.05, 0.99)
    relative_floor: float = 1.0e-12
    relative_valence_weight: float = 0.0
    relative_valence_x_range: tuple[float, float] = (1.0e-4, 0.1)
    relative_valence_floor: float = 1.0e-12
    relative_f2_weight: float = 0.0
    relative_f2_x_range: tuple[float, float] = (1.0e-4, 0.1)
    relative_f2_floor: float = 1.0e-12

    @classmethod
    def build(
        cls,
        basis_set: str,
        n_basis: int,
        flavors: Sequence[int],
        x_grid: Sequence[float],
        q: float,
        metric: str,
        **kwargs: object,
    ) -> "Figure2ProjectionOperator":
        if metric != "relative_gluon":
            return cls(ProjectionOperator.build(basis_set, n_basis, flavors, x_grid, q, metric), metric)
        if 21 not in flavors:
            raise ValueError("relative_gluon requires PID 21 in flavors")
        values = {
            name: kwargs.get(name, default)
            for name, default in (
                ("relative_weight", 0.1),
                ("relative_x_range", (0.05, 0.99)),
                ("relative_floor", 1.0e-12),
                ("relative_valence_weight", 0.0),
                ("relative_valence_x_range", (1.0e-4, 0.1)),
                ("relative_valence_floor", 1.0e-12),
                ("relative_f2_weight", 0.0),
                ("relative_f2_x_range", (1.0e-4, 0.1)),
                ("relative_f2_floor", 1.0e-12),
            )
        }
        for name in ("relative_weight", "relative_floor", "relative_valence_floor", "relative_f2_floor"):
            if float(values[name]) <= 0:
                raise ValueError(f"{name} must be positive")
        for name in ("relative_x_range", "relative_valence_x_range", "relative_f2_x_range"):
            low, high = (float(v) for v in values[name])
            if not 0 < low < high <= 1:
                raise ValueError(f"{name} must satisfy 0 < low < high <= 1")
            values[name] = (low, high)
        for name in ("relative_valence_weight", "relative_f2_weight"):
            if float(values[name]) < 0:
                raise ValueError(f"{name} must be non-negative")
        return cls(
            ProjectionOperator.build(basis_set, n_basis, flavors, x_grid, q, "dist0"),
            metric,
            **{name: float(value) if "range" not in name else value for name, value in values.items()},
        )

    def __getattr__(self, name: str) -> object:
        return getattr(self.base, name)

    def project_grid(self, target_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.metric != "relative_gluon":
            return self.base.project_grid(target_grid)
        target_grid = np.asarray(target_grid, dtype=float)
        if target_grid.shape != self.base.reference_grid.shape:
            raise ValueError(f"Expected target grid shape {self.base.reference_grid.shape}, got {target_grid.shape}")
        displacement = (target_grid - self.base.reference_grid).reshape(-1)
        matrix_3d = self.base.matrix.reshape(len(self.base.basis.flavors), len(self.base.x_grid), -1)
        displacement_2d = displacement.reshape(len(self.base.basis.flavors), len(self.base.x_grid))
        matrices = [self.base.matrix]
        displacements = [displacement]

        def add_relative(matrix: np.ndarray, target: np.ndarray, displacement: np.ndarray, mask: np.ndarray, weight: float, floor: float) -> None:
            denominator = np.maximum(np.abs(target[mask]), floor)
            matrices.append(weight * matrix[mask] / denominator[:, None])
            displacements.append(weight * displacement[mask] / denominator)

        gluon = self.base.basis.flavors.index(21)
        gluon_mask = (self.base.x_grid >= self.relative_x_range[0]) & (self.base.x_grid <= self.relative_x_range[1])
        add_relative(matrix_3d[gluon], target_grid[gluon], displacement_2d[gluon], gluon_mask, self.relative_weight, self.relative_floor)

        if self.relative_valence_weight:
            mask = (self.base.x_grid >= self.relative_valence_x_range[0]) & (self.base.x_grid <= self.relative_valence_x_range[1])
            for quark, antiquark in ((2, -2), (1, -1)):
                q, aq = self.base.basis.flavors.index(quark), self.base.basis.flavors.index(antiquark)
                add_relative(matrix_3d[q] - matrix_3d[aq], target_grid[q] - target_grid[aq], displacement_2d[q] - displacement_2d[aq], mask, self.relative_valence_weight, self.relative_valence_floor)

        if self.relative_f2_weight:
            charges = {2: 4 / 9, -2: 4 / 9, 1: 1 / 9, -1: 1 / 9, 3: 1 / 9, -3: 1 / 9, 4: 4 / 9, -4: 4 / 9, 5: 1 / 9, -5: 1 / 9}
            missing = set(charges) - set(self.base.basis.flavors)
            if missing:
                raise ValueError(f"relative_f2 requires flavors {sorted(missing)}")
            mask = (self.base.x_grid >= self.relative_f2_x_range[0]) & (self.base.x_grid <= self.relative_f2_x_range[1])
            def weighted(values: np.ndarray) -> np.ndarray:
                return sum(charge * values[self.base.basis.flavors.index(pid)] for pid, charge in charges.items())
            add_relative(weighted(matrix_3d), weighted(target_grid), weighted(displacement_2d), mask, self.relative_f2_weight, self.relative_f2_floor)

        coefficients, *_ = np.linalg.lstsq(np.vstack(matrices), np.concatenate(displacements), rcond=None)
        projected = self.base.reference_grid + (self.base.matrix @ coefficients).reshape(self.base.reference_grid.shape)
        return projected, coefficients, target_grid - projected
