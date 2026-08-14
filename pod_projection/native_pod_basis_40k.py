#!/usr/bin/env python3
"""Load native shifts from the private 250503 POD LHAPDF set.

The basis vectors are used without additional normalization:

    basis_i(x, Q, pid) = member_i(x, Q, pid) - member_0(x, Q, pid)

This is the loader used by the GOLLUM projection demo, kept standalone so it
can also be reused by a future xFitter interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import lhapdf
import numpy as np


PDF_SET = "250503_pod_basis_40k"


def pdfset_size(pdf_set: str) -> int:
    """Return the number of members for an installed LHAPDF set."""
    size = lhapdf.getPDFSet(pdf_set).size
    return int(size() if callable(size) else size)


@dataclass(frozen=True)
class NativePODBasis:
    pdf_set: str
    reference_pdf: object
    basis_pdfs: Sequence[object]
    variations: Sequence[int]
    flavors: Sequence[int]

    @classmethod
    def load(
        cls,
        pdf_set: str = PDF_SET,
        variations: Iterable[int] | None = None,
        flavors: Iterable[int] | None = None,
    ) -> "NativePODBasis":
        n_members = pdfset_size(pdf_set)

        if variations is None:
            variations = range(1, n_members)
        variations = tuple(int(member) for member in variations)
        invalid = [member for member in variations if not 1 <= member < n_members]
        if invalid:
            raise ValueError(
                f"Basis members must lie in 1..{n_members - 1}; got {invalid}"
            )

        reference_pdf = lhapdf.mkPDF(pdf_set, 0)
        basis_pdfs = tuple(lhapdf.mkPDF(pdf_set, member) for member in variations)

        if flavors is None:
            flavors = reference_pdf.flavors()
        flavors = tuple(int(pid) for pid in flavors)

        return cls(
            pdf_set=pdf_set,
            reference_pdf=reference_pdf,
            basis_pdfs=basis_pdfs,
            variations=variations,
            flavors=flavors,
        )

    @property
    def nvariations(self) -> int:
        return len(self.variations)

    def xfx_grid(self, pdf: object, x_grid: Sequence[float], q: float) -> np.ndarray:
        """Evaluate x*f(x,Q), returning an array with shape (flavor, x)."""
        x_grid = np.asarray(x_grid, dtype=float)
        values = np.empty((len(self.flavors), len(x_grid)), dtype=float)
        for i_pid, pid in enumerate(self.flavors):
            values[i_pid] = [pdf.xfxQ(pid, x, q) for x in x_grid]
        return values

    def reference_grid(self, x_grid: Sequence[float], q: float) -> np.ndarray:
        return self.xfx_grid(self.reference_pdf, x_grid, q)

    def native_shift_grid(self, x_grid: Sequence[float], q: float) -> np.ndarray:
        """Return native member_i-member_0 shifts as (mode, flavor, x)."""
        reference = self.reference_grid(x_grid, q)
        shifts = np.empty(
            (self.nvariations, len(self.flavors), len(x_grid)), dtype=float
        )
        for i_member, pdf in enumerate(self.basis_pdfs):
            shifts[i_member] = self.xfx_grid(pdf, x_grid, q) - reference
        return shifts

    def combine(
        self, coeffs: Sequence[float], x_grid: Sequence[float], q: float
    ) -> np.ndarray:
        """Construct member_0 + sum_i(coeff_i * basis_i)."""
        coeffs = np.asarray(coeffs, dtype=float)
        if coeffs.shape != (self.nvariations,):
            raise ValueError(
                f"Expected {self.nvariations} coefficients, got {coeffs.shape}"
            )
        reference = self.reference_grid(x_grid, q)
        shifts = self.native_shift_grid(x_grid, q)
        return reference + np.einsum("i,ifx->fx", coeffs, shifts)


if __name__ == "__main__":
    basis = NativePODBasis.load()
    print(f"Loaded {basis.pdf_set}: {basis.nvariations} native shifts")
    print(f"Flavors: {list(basis.flavors)}")
