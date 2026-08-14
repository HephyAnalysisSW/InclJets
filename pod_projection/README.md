# External PDF projection onto the native POD basis

This directory is a self-contained numerical version of the GOLLUM CT18NNLO
projection demo. It uses LHAPDF, NumPy, and Matplotlib; it does not require the
site-specific GOLLUM `common` modules or its plot/synchronization setup.

The implementation follows
[`plot_ct18_pod_projection_uncertainties.py`](https://github.com/HephyAnalysisSW/GOLLUM/blob/746fc21ff8fd223e656bfe979813c061657805c4/plot/pod/plot_ct18_pod_projection_uncertainties.py),
[`check_native_pod_completeness.py`](https://github.com/HephyAnalysisSW/GOLLUM/blob/746fc21ff8fd223e656bfe979813c061657805c4/plot/pod/check_native_pod_completeness.py),
and
[`native_pod_basis_40k.py`](https://github.com/HephyAnalysisSW/GOLLUM/blob/746fc21ff8fd223e656bfe979813c061657805c4/plot/pod/native_pod_basis_40k.py)
at GOLLUM commit `746fc21ff8fd223e656bfe979813c061657805c4`.

## Quick start

Activate the Miniforge environment that contains LHAPDF, then run:

```bash
cd pod_projection
./run_ct18_demo.sh
```

The first command inside the wrapper projects all 59 CT18NNLO members and
writes `outputs/ct18nnlo_pod_projection.npz`. The second makes PDF and PNG
comparisons under `plots/`. The `ratio_bands` output reproduces the upstream
observable: the CT18 and projected central values and 90% bands are all divided
by the CT18 central PDF. The separate `ratios` output is a reconstruction
diagnostic comparing projected/target central values and uncertainties.

To project another installed LHAPDF set:

```bash
python project_pdf.py --target-set PDF4LHC21_40_pdfas \
  --output outputs/pdf4lhc21_pod_projection.npz
```

For a different basis size, flavor selection, integration metric, or grid
slice, use `python project_pdf.py --help`. Named flavor selections are `qcd5`,
`qcd4`, and `paper`; a comma-separated PID list is also accepted.

## Numerical definition

The reference is member 0 of `250503_pod_basis_40k`. Mode `i` is the native
LHAPDF shift `member_i - member_0`, without a `max_amplitudes` rescaling. At
the defaults, each target member is sampled as `x*f(x,Q)` at `Q=1.65 GeV` for
11 QCD5 flavors and 140 x points (`LHAPDF_XGRID[36:-20]`). The 100 coefficients
solve

```text
(X^T W X) c = X^T W (target - reference),
projected = reference + X c.
```

The upstream default metric is `dist0`, so `W` is the identity. For a paired
Hessian target such as CT18NNLO, the native symmetric uncertainty is
`0.5*sqrt(sum_k(member_k+ - member_k-)^2)`. CT18NNLO supplies 29 eigenvector
pairs at 90% confidence level; no 68% rescaling is applied here.

## Saved arrays

The NPZ file contains the x grid, flavors, POD reference and shifts, external
and reconstructed member grids, residuals, one coefficient vector per member,
Hessian bands, coefficient covariance/displacements, the Gram matrix, and JSON
metadata. This is the natural intermediate format for the later xFitter
interface; it keeps the projection definition separate from xFitter's fit and
likelihood machinery.

## Current CT18 validation

The checked 100-mode run has 1540 rows, rank 100, and reconstructs the CT18NNLO
central displacement with a weighted residual/shift ratio of `5.4828e-4`. The
median across all 59 members is `6.310e-4`. The 29-pair coefficient covariance
has rank 29, as it must.

The native coordinates are numerically anisotropic: `cond(X)=2.20e7` and
`cond(X^T X)=4.86e14`. This code deliberately retains GOLLUM's normal-equation
solve. A separate SVD least-squares audit changed the CT18 central coefficients
by `3.4e-8` relatively and the reconstructed shift by only `9.4e-12`, so the
present projection is stable. The conditioning should nevertheless be handled
explicitly when the POD coefficients become xFitter parameters.
