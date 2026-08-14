#!/usr/bin/env python3
"""Numerical projection of an external LHAPDF set onto a native POD basis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import lhapdf
import numpy as np

try:
    from .native_pod_basis_40k import NativePODBasis, pdfset_size
except ImportError:  # Support direct execution from this directory.
    from native_pod_basis_40k import NativePODBasis, pdfset_size


Q0 = 1.65
# Flavor order used by the upstream plotting demo (not the differently ordered
# tuple in its numerical completeness checker).
QCD5_FLAVORS = (21, 2, -2, 1, -1, 3, -3, 4, -4, 5, -5)
QCD4_FLAVORS = (21, 2, -2, 1, -1, 3, -3, 4, -4)
PAPER_FLAVORS = (1, -1, 2, -2, 3, -3, 4, 21)
METRICS = ("dist0", "dist4_x", "trapz_x", "trapz_logx", "relative_gluon")

# Exact grid used by the upstream GOLLUM demo.  It is the 196-point x grid of
# the private POD LHAPDF set, retained here at the precision of that demo.
LHAPDF_XGRID = np.asarray([1e-09, 1.29708482343957e-09, 1.68242903474257e-09, 2.18225315420583e-09, 2.83056741739819e-09, 3.67148597892941e-09, 4.76222862935315e-09, 6.1770142737618e-09, 8.01211109898438e-09, 1.03923870607245e-08, 1.34798064073805e-08, 1.74844503691778e-08, 2.26788118881103e-08, 2.94163370300835e-08, 3.81554746595878e-08, 4.94908707232129e-08, 6.41938295708371e-08, 8.32647951986859e-08, 1.08001422993829e-07, 1.4008687308113e-07, 1.81704331793772e-07, 2.35685551545377e-07, 3.05703512595323e-07, 3.96522309841747e-07, 5.1432125723657e-07, 6.67115245136676e-07, 8.65299922973143e-07, 1.12235875241487e-06, 1.45577995547683e-06, 1.88824560514613e-06, 2.44917352454946e-06, 3.17671650028717e-06, 4.12035415232797e-06, 5.3442526575209e-06, 6.93161897806315e-06, 8.99034258238145e-06, 1.16603030112258e-05, 1.51228312288769e-05, 1.96129529349212e-05, 2.54352207134502e-05, 3.29841683435992e-05, 4.27707053972016e-05, 5.54561248105849e-05, 7.18958313632514e-05, 9.31954227979614e-05, 0.00012078236773133, 0.000156497209466554, 0.000202708936328495, 0.000262459799331951, 0.000339645244168985, 0.000439234443000422, 0.000567535660104533, 0.000732507615725537, 0.000944112105452451, 0.00121469317686978, 0.00155935306118224, 0.00199627451141338, 0.00254691493736552, 0.00323597510213126, 0.00409103436509565, 0.00514175977083962, 0.00641865096062317, 0.00795137940306351, 0.009766899996241, 0.0118876139251364, 0.0143298947643919, 0.0171032279460271, 0.0202100733925079, 0.0236463971369542, 0.0274026915728357, 0.0314652506132444, 0.0358174829282429, 0.0404411060163317, 0.0453171343973807, 0.0504266347950069, 0.0557512610084339, 0.0612736019390519, 0.0669773829498255, 0.0728475589986517, 0.0788703322292727, 0.0850331197801452, 0.0913244910278679, 0.0977340879783772, 0.104252538208639, 0.110871366547237, 0.117582909372878, 0.124380233801599, 0.131257062945031, 0.138207707707289, 0.145227005135651, 0.152310263065985, 0.159453210652156, 0.166651954293987, 0.173902938455578, 0.181202910873333, 0.188548891679097, 0.195938145999193, 0.203368159629765, 0.210836617429103, 0.218341384106561, 0.225880487124065, 0.233452101459503, 0.241054536011681, 0.248686221452762, 0.256345699358723, 0.264031612468684, 0.271742695942783, 0.279477769504149, 0.287235730364833, 0.295015546847664, 0.302816252626866, 0.310636941519503, 0.318476762768082, 0.326334916761672, 0.334210651149156, 0.342103257303627, 0.350012067101685, 0.357936449985571, 0.365875810279643, 0.373829584735962, 0.381797240286494, 0.389778271981947, 0.397772201099286, 0.40577857340234, 0.413796957540671, 0.421826943574548, 0.429868141614175, 0.437920180563205, 0.44598270695699, 0.454055383887562, 0.462137890007651, 0.470229918607142, 0.478331176755675, 0.486441384506059, 0.494560274153348, 0.502687589545177, 0.510823085439086, 0.518966526903235, 0.527117688756998, 0.535276355048428, 0.543442318565661, 0.551615380379768, 0.559795349416641, 0.5679820420558, 0.576175281754088, 0.584374898692498, 0.59258072944444, 0.60079261666395, 0.609010408792398, 0.61723395978245, 0.625463128838069, 0.633697780169485, 0.641937782762089, 0.650183010158361, 0.658433340251944, 0.666688655093089, 0.674948840704708, 0.683213786908386, 0.691483387159697, 0.699757538392251, 0.708036140869916, 0.716319098046733, 0.724606316434025, 0.732897705474271, 0.741193177421404, 0.749492647227008, 0.757796032432224, 0.766103253064927, 0.774414231541921, 0.782728892575836, 0.791047163086478, 0.799368972116378, 0.807694250750291, 0.816022932038457, 0.824354950923382, 0.832690244169987, 0.841028750298844, 0.8493704095226, 0.857715163684985, 0.866062956202683, 0.874413732009721, 0.882767437504206, 0.891124020497459, 0.899483430165226, 0.907845617001021, 0.916210532771399, 0.924578130473112, 0.932948364292029, 0.941321189563734, 0.949696562735755, 0.958074441331298, 0.966454783914439, 0.974837550056705, 0.983222700304978, 0.991610196150662, 1.0], dtype=float)


def parse_flavors(value: str, basis_set: str, target_set: str) -> tuple[int, ...]:
    """Parse a named flavor selection or a comma-separated list of PIDs."""
    if value == "qcd5":
        requested = QCD5_FLAVORS
    elif value == "qcd4":
        requested = QCD4_FLAVORS
    elif value == "paper":
        requested = PAPER_FLAVORS
    else:
        requested = tuple(int(pid.strip()) for pid in value.split(",") if pid.strip())
    if not requested:
        raise ValueError("No flavors selected")

    basis_flavors = set(lhapdf.mkPDF(basis_set, 0).flavors())
    target_flavors = set(lhapdf.mkPDF(target_set, 0).flavors())
    missing = [
        pid for pid in requested if pid not in basis_flavors or pid not in target_flavors
    ]
    if missing:
        raise ValueError(f"Requested PIDs unavailable in both PDF sets: {missing}")
    return tuple(requested)


def trapezoid_weights(grid: Sequence[float]) -> np.ndarray:
    grid = np.asarray(grid, dtype=float)
    if grid.ndim != 1 or len(grid) < 2:
        raise ValueError("Expected a one-dimensional grid with at least two points")
    weights = np.empty_like(grid)
    weights[0] = 0.5 * (grid[1] - grid[0])
    weights[-1] = 0.5 * (grid[-1] - grid[-2])
    weights[1:-1] = 0.5 * (grid[2:] - grid[:-2])
    return weights


def metric_weights(metric: str, x_grid: np.ndarray, n_flavors: int) -> np.ndarray:
    """Build the diagonal metric used in the least-squares projection."""
    if metric == "dist0":
        return np.ones(n_flavors * len(x_grid))
    if metric == "dist4_x":
        return np.tile(np.abs(x_grid), n_flavors)
    if metric == "trapz_x":
        return np.tile(trapezoid_weights(x_grid), n_flavors)
    if metric == "trapz_logx":
        return np.tile(trapezoid_weights(np.log(x_grid)), n_flavors)
    raise ValueError(f"Unknown metric {metric!r}; choose one of {METRICS}")


@dataclass(frozen=True)
class ProjectionOperator:
    """Precomputed weighted least-squares operator for one POD configuration."""

    basis: NativePODBasis
    x_grid: np.ndarray
    q: float
    metric: str
    reference_grid: np.ndarray
    shift_grid: np.ndarray
    matrix: np.ndarray
    weights: np.ndarray
    gram: np.ndarray
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
        q: float = Q0,
        metric: str = "dist0",
        relative_weight: float = 0.1,
        relative_x_range: tuple[float, float] = (0.05, 0.99),
        relative_floor: float = 1.0e-12,
        relative_valence_weight: float = 0.0,
        relative_valence_x_range: tuple[float, float] = (1.0e-4, 0.1),
        relative_valence_floor: float = 1.0e-12,
        relative_f2_weight: float = 0.0,
        relative_f2_x_range: tuple[float, float] = (1.0e-4, 0.1),
        relative_f2_floor: float = 1.0e-12,
    ) -> "ProjectionOperator":
        if n_basis < 1 or n_basis >= pdfset_size(basis_set):
            raise ValueError(
                f"n_basis must be in 1..{pdfset_size(basis_set) - 1}"
            )
        x_grid = np.asarray(x_grid, dtype=float)
        basis = NativePODBasis.load(
            basis_set,
            variations=range(1, n_basis + 1),
            flavors=flavors,
        )
        reference = basis.reference_grid(x_grid, q)
        shifts = basis.native_shift_grid(x_grid, q)
        matrix = shifts.reshape(n_basis, -1).T
        if metric == "relative_gluon":
            # The absolute part remains unweighted; the gluon-relative rows
            # are assembled in project_grid because their denominator is the
            # external target PDF and therefore varies point by point.
            weights = metric_weights("dist0", x_grid, len(flavors))
            if 21 not in flavors:
                raise ValueError("relative_gluon requires PID 21 in flavors")
            if relative_weight <= 0:
                raise ValueError("relative_weight must be positive")
            if not (0 < relative_x_range[0] < relative_x_range[1] <= 1):
                raise ValueError("relative_x_range must satisfy 0 < low < high <= 1")
            if relative_floor <= 0:
                raise ValueError("relative_floor must be positive")
            if relative_valence_weight < 0:
                raise ValueError("relative_valence_weight must be non-negative")
            if not (
                0 < relative_valence_x_range[0]
                < relative_valence_x_range[1]
                <= 1
            ):
                raise ValueError(
                    "relative_valence_x_range must satisfy 0 < low < high <= 1"
                )
            if relative_valence_floor <= 0:
                raise ValueError("relative_valence_floor must be positive")
            if relative_f2_weight < 0:
                raise ValueError("relative_f2_weight must be non-negative")
            if not (0 < relative_f2_x_range[0] < relative_f2_x_range[1] <= 1):
                raise ValueError("relative_f2_x_range must satisfy 0 < low < high <= 1")
            if relative_f2_floor <= 0:
                raise ValueError("relative_f2_floor must be positive")
        else:
            weights = metric_weights(metric, x_grid, len(flavors))
        gram = matrix.T @ (weights[:, np.newaxis] * matrix)
        return cls(
            basis=basis,
            x_grid=x_grid,
            q=float(q),
            metric=metric,
            reference_grid=reference,
            shift_grid=shifts,
            matrix=matrix,
            weights=weights,
            gram=gram,
            relative_weight=float(relative_weight),
            relative_x_range=tuple(float(v) for v in relative_x_range),
            relative_floor=float(relative_floor),
            relative_valence_weight=float(relative_valence_weight),
            relative_valence_x_range=tuple(
                float(v) for v in relative_valence_x_range
            ),
            relative_valence_floor=float(relative_valence_floor),
            relative_f2_weight=float(relative_f2_weight),
            relative_f2_x_range=tuple(float(v) for v in relative_f2_x_range),
            relative_f2_floor=float(relative_f2_floor),
        )

    @property
    def rank(self) -> int:
        return int(np.linalg.matrix_rank(self.matrix))

    @property
    def condition_number(self) -> float:
        return float(np.linalg.cond(self.gram))

    def evaluate(self, pdf: object) -> np.ndarray:
        return self.basis.xfx_grid(pdf, self.x_grid, self.q)

    def project_grid(
        self, target_grid: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return projected grid, POD coefficients, and target-project residual."""
        target_grid = np.asarray(target_grid, dtype=float)
        if target_grid.shape != self.reference_grid.shape:
            raise ValueError(
                f"Expected target grid shape {self.reference_grid.shape}, "
                f"got {target_grid.shape}"
            )
        displacement = (target_grid - self.reference_grid).reshape(-1)
        if self.metric == "relative_gluon":
            # Add lambda * (target-project)/target for the selected gluon
            # x-range. This is a configurable closure metric, not an
            # experimental likelihood term or a PDF prior.
            augmented_matrix = [self.matrix]
            augmented_displacement = [displacement]
            gluon = self.basis.flavors.index(21)
            target_gluon = target_grid[gluon]
            mask = (self.x_grid >= self.relative_x_range[0]) & (
                self.x_grid <= self.relative_x_range[1]
            )
            denominator = np.maximum(np.abs(target_gluon[mask]), self.relative_floor)
            relative_matrix = self.matrix.reshape(len(self.basis.flavors), len(self.x_grid), -1)[gluon, mask]
            relative_displacement = displacement.reshape(len(self.basis.flavors), len(self.x_grid))[gluon, mask]
            augmented_matrix.append(self.relative_weight * relative_matrix / denominator[:, None])
            augmented_displacement.append(self.relative_weight * relative_displacement / denominator)
            if self.relative_valence_weight:
                valence_mask = (
                    (self.x_grid >= self.relative_valence_x_range[0])
                    & (self.x_grid <= self.relative_valence_x_range[1])
                )
                for quark, antiquark in ((2, -2), (1, -1)):
                    q_index = self.basis.flavors.index(quark)
                    aq_index = self.basis.flavors.index(antiquark)
                    valence_target = target_grid[q_index, valence_mask] - target_grid[
                        aq_index, valence_mask
                    ]
                    valence_matrix = (
                        self.matrix.reshape(
                            len(self.basis.flavors), len(self.x_grid), -1
                        )[q_index, valence_mask]
                        - self.matrix.reshape(
                            len(self.basis.flavors), len(self.x_grid), -1
                        )[aq_index, valence_mask]
                    )
                    valence_displacement = displacement.reshape(
                        len(self.basis.flavors), len(self.x_grid)
                    )[q_index, valence_mask] - displacement.reshape(
                        len(self.basis.flavors), len(self.x_grid)
                    )[aq_index, valence_mask]
                    valence_denominator = np.maximum(
                        np.abs(valence_target), self.relative_valence_floor
                    )
                    augmented_matrix.append(
                        self.relative_valence_weight
                        * valence_matrix
                        / valence_denominator[:, None]
                    )
                    augmented_displacement.append(
                        self.relative_valence_weight
                        * valence_displacement
                        / valence_denominator
                    )
            if self.relative_f2_weight:
                # Photon-exchange F2 proxy: sum_f e_f^2 x(q_f + qbar_f).
                # It is a projection constraint for the HERA-sensitive light
                # sea combination, not an extra term in the fit likelihood.
                f2_mask = (
                    (self.x_grid >= self.relative_f2_x_range[0])
                    & (self.x_grid <= self.relative_f2_x_range[1])
                )
                f2_components = {
                    2: 4.0 / 9.0,
                    -2: 4.0 / 9.0,
                    1: 1.0 / 9.0,
                    -1: 1.0 / 9.0,
                    3: 1.0 / 9.0,
                    -3: 1.0 / 9.0,
                    4: 4.0 / 9.0,
                    -4: 4.0 / 9.0,
                    5: 1.0 / 9.0,
                    -5: 1.0 / 9.0,
                }
                missing = set(f2_components) - set(self.basis.flavors)
                if missing:
                    raise ValueError(f"relative_f2 requires flavors {sorted(missing)}")
                matrix_3d = self.matrix.reshape(
                    len(self.basis.flavors), len(self.x_grid), -1
                )
                displacement_2d = displacement.reshape(
                    len(self.basis.flavors), len(self.x_grid)
                )
                f2_matrix = sum(
                    charge * matrix_3d[self.basis.flavors.index(pid), f2_mask]
                    for pid, charge in f2_components.items()
                )
                f2_displacement = sum(
                    charge * displacement_2d[self.basis.flavors.index(pid), f2_mask]
                    for pid, charge in f2_components.items()
                )
                f2_target = sum(
                    charge * target_grid[self.basis.flavors.index(pid), f2_mask]
                    for pid, charge in f2_components.items()
                )
                f2_denominator = np.maximum(
                    np.abs(f2_target), self.relative_f2_floor
                )
                augmented_matrix.append(
                    self.relative_f2_weight * f2_matrix / f2_denominator[:, None]
                )
                augmented_displacement.append(
                    self.relative_f2_weight * f2_displacement / f2_denominator
                )
            coefficients, *_ = np.linalg.lstsq(
                np.vstack(augmented_matrix), np.concatenate(augmented_displacement), rcond=None
            )
        else:
            rhs = self.matrix.T @ (self.weights * displacement)
            # Keep the same normal-equation solve as the upstream demonstration.
            coefficients = np.linalg.solve(self.gram, rhs)
        projected = self.reference_grid + (
            self.matrix @ coefficients
        ).reshape(self.reference_grid.shape)
        return projected, coefficients, target_grid - projected


def hessian_symmetric_band(member_grids: np.ndarray) -> np.ndarray:
    """Symmetric Hessian band: 1/2 sqrt(sum_k(PDF_k+ - PDF_k-)^2)."""
    member_grids = np.asarray(member_grids, dtype=float)
    if member_grids.shape[0] < 3 or (member_grids.shape[0] - 1) % 2:
        raise ValueError("Expected central plus an even number of Hessian members")
    differences = member_grids[1::2] - member_grids[2::2]
    return 0.5 * np.sqrt(np.sum(differences * differences, axis=0))


def coefficient_covariance(coefficients: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hessian coefficient covariance and pairwise displacements."""
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape[0] < 3 or (coefficients.shape[0] - 1) % 2:
        raise ValueError("Expected central plus an even number of Hessian members")
    displacements = 0.5 * (coefficients[1::2] - coefficients[2::2])
    return displacements.T @ displacements, displacements
