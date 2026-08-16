# CMS-only POD closure diagnostic

This diagnostic keeps the direct Figure-2 PDF point and its full-100-mode POD
projection fixed. It removes HERA completely, selects the four CMS inclusive-
jet rapidity datasets, and fixes all 29 active CMS-only correlated-nuisance
shifts to zero. (The combined reference bookkeeping has 36 CMS-labelled
entries; the seven combined-fit `proc_*` terms are not in the CMS-only active
Hessian source list.)
There is no minimization or nuisance profiling.

After a completed full-POD likelihood scan, run:

```bash
cd fits/SMP-22-014/figure2/cms_closure
python run_zero_systematics.py \
  --source ../likelihood_scans/production_16x21_gluon_valence_f2 \
  --output output
```

For the PDF plot, use the full-POD projection from the selected scan reference,
for example `../likelihood_scans/production_16x21_gluon_valence_f2/runs/_reference/full_pod_projection.npz`:

```bash
python plot_closure.py \
  --projection ../likelihood_scans/production_16x21_gluon_valence_f2/runs/_reference/full_pod_projection.npz \
  --cms output/cms_zero_systematics.npz --output-dir output/plots
```

The PDF figure uses the signed POD-direct difference normalised to the maximum
absolute direct PDF of each flavour. This avoids undefined ratios at flavour
zero crossings near \(x=1\) while ranking the largest meaningful defect by
flavour and \(x\). The prediction figure shows data, direct theory, POD theory,
and their bin-by-bin ratio. The residual figure is a diagonal uncorrelated-residual
*proxy*, useful to localise bins but not a replacement for xFitter's exact
Poisson/Hessian likelihood (reported in `metadata.yaml`). The experimental
correlated systematic-nuisance shifts are fixed to zero.
statistical covariance supplied with the CMS data remains active; only the 29
active correlated systematic-nuisance shifts are fixed to zero.

## Gluon-only isolation

The following control constructs a full direct input grid down to \(x=10^{-9}\),
then evaluates two otherwise identical QCDNUM evolutions at \(Q_0=1.65\) GeV:
the direct grid and a hybrid in which only the gluon is replaced by its full-POD
reconstruction. It is therefore insensitive to residual POD differences in the
quark flavours.

```bash
python run_gluon_isolation.py \
  --source ../likelihood_scans/production_16x21_gluon_valence_f2 \
  --output output_gluon_isolation
```

Read `output_gluon_isolation/summary.yaml`; the field
`total_chi2_shifts.gluon_hybrid` is the desired result. The direct-input control
is deliberately not compared directly with the native HERAPDF evaluation,
since it starts evolution at 1.65 GeV rather than the native 1.378 GeV. It is
compared only with the hybrid, which shares exactly the same LHAPDF grid and
evolution route.

To replace each flavour separately and test the additive accounting, add
`--all-flavors`. The summary contains `total_chi2_shifts`, including the
all-flavour hybrid. Make the attribution figure with:

```bash
python run_gluon_isolation.py --all-flavors \
  --source ../likelihood_scans/production_16x21_gluon_valence_f2 \
  --output output_flavour_isolation
python plot_flavour_attribution.py output_flavour_isolation/summary.yaml \
  --output output_flavour_isolation/flavour_attribution.png
```

For the nominal POD-closure comparison, use the faster scale-aligned mode:

```bash
python run_gluon_isolation.py --all-flavours-only \
  --source ../likelihood_scans/production_16x21_gluon_valence_f2 \
  --output output_scale_aligned_closure
```

Its `scale_aligned_full_POD_closure_delta_chi2` compares the direct and POD
PDFs with the same full LHAPDF grid, the same \(Q_0=1.65\) GeV, and the same
QCDNUM evolution. This is the POD-closure result; do not use the native
HERAPDF-to-POD difference for that purpose.

Make the three publication-oriented closure figures after running both the
scale-aligned and technical controls:

```bash
python plot_scale_aligned_closure.py \
  --projection ../likelihood_scans/production_16x21_gluon_valence_f2/runs/_reference/full_pod_projection.npz \
  --closure output_scale_aligned_closure/summary.yaml \
  --closure-data output_scale_aligned_closure/flavour_isolation.npz \
  --technical output_technical_controls/summary.yaml \
  --output-dir output_scale_aligned_closure/plots
```

## Technical controls: grid versus input scale

The matched flavour study intentionally starts QCDNUM at 1.65 GeV, whereas
the native HERAPDF parametrisation starts at 1.378404875209 GeV. Separate the
pure LHAPDF-table import from this change of input scale with:

```bash
python run_technical_controls.py \
  --source ../likelihood_scans/production_16x21_gluon_valence_f2 \
  --output output_technical_controls
```

It evaluates native HERAPDF, the same full \(x\)-grid table re-imported at its
native input scale, and the table re-imported at 1.65 GeV. `summary.yaml`
reports `native_to_table_at_native_q0` (grid/interpolation) separately from
`table_native_q0_to_table_q165` (input-scale/evolution-history change).
correlated systematic-nuisance shifts are fixed to zero.
