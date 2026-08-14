# Seven-mode gluon-POD CMS test

This is a prepared, deliberately unrun CMS-inclusive-jet-only fit using modes
1--7 of `gluon_POD_nongluon_PDF4LHC21`. Its `pod.yaml` sets both the external
projection and the xFitter input deformation to PID 21; the decomposition
therefore retains the reference member's non-gluon flavors exactly at Q0.

Run the PDF-level basis closure test with:

```bash
./validate_basis.sh
```

The validation checks the native member convention and loads this exact card
through xFitter, then checks linear QCDNUM evolution at Q = 1.65, 10, 100, and
1000 GeV. No data, likelihood, or minimization is invoked. `run.sh` requires
the explicit `--run-fit` argument after the cards have been reviewed.

The stored LHAPDF modes are gluon-only up to table precision. Modes 4--7 have
non-gluon residues no larger than 7.2e-7 in x*f on the projection grid (below
9.2e-7 of their gluon shift). `active_flavors: [21]` removes even that residue
from the xFitter input definition. Change `basis_members` to alter the fitted
truncation; omit it or use `all` to select every non-reference member.
