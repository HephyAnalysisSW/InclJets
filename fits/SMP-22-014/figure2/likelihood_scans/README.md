# Direct versus full-POD likelihood scans

This directory evaluates the Figure 2 likelihood without a fit. At every scan
point all 16 external HERAPDF parameters are fixed, all 198 nuisance shifts are
fixed to the stored global-fit values, and xFitter is called once (`call fcn 3`).
There is no minimization and no profiling.

The comparison follows this chain:

1. xFitter evaluates the direct HERAPDF likelihood and exports the evolved PDF
   at `Q_ext = 1.65 GeV` on the configured POD projection x grid (currently
   the 160 points from `x=0.0001208` through `x=1`).
2. Python projects all 11 QCD5 flavors on all 100 native directions of
   `250503_pod_basis_40k` using the configured metric. The current
   `relative_gluon` metric augments the absolute residual with
   `lambda*(target gluon - POD gluon)/target gluon` for `0.05<=x<=0.99`.
3. xFitter evaluates the same fixed-nuisance likelihood with those 100 POD
   coefficients and the scan point's separately supplied alpha_s.
4. The runner stores compact row arrays in `scan_results.npz`; plotting reads
   only this file and never invokes xFitter.

The HERA and CMS terms each include their data, log-penalty, and assigned
correlated-nuisance penalty contributions. Their sum exactly reproduces the
joint-card total.

## Quick validated example

The checked three-point `B_g` scan and plots are in `smoke_Bg_3point`:

```bash
cd fits/SMP-22-014/figure2/likelihood_scans
python run_likelihood_scan.py \
  --output smoke_Bg_3point \
  --parameters Bg \
  --coordinates=-0.5,0,0.5
python plot_likelihood_scan.py \
  --input smoke_Bg_3point/scan_results.npz \
  --output-dir smoke_Bg_3point/plots
```

Existing completed point directories are reused, so the runner can be invoked
again after an interruption. A point is never overwritten by the underlying
evaluators. Incomplete point directories from an interrupted evaluator are
removed automatically; completed directories remain protected. Do not run two
instances against the same output directory concurrently.

## Full configured scan

Omitting the parameter and coordinate selectors uses all 16 independent
parameters and the 21 points from `scan_config.yaml`:

```bash
cd fits/SMP-22-014/figure2/likelihood_scans
python run_likelihood_scan.py --output production_16x21_relative_gluon
python plot_likelihood_scan.py \
  --input production_16x21_relative_gluon/scan_results.npz \
  --output-dir production_16x21_relative_gluon/plots
```

On this machine a direct point takes about 13 seconds and a full-POD point
about 17 seconds. Reusing the common central point, the full 16x21 scan is
expected to take roughly 2.8 hours serially and about 0.4 GB. The previous
`production_16x21` run used the old 140-node `dist0` metric and was stopped;
its results remain as a historical diagnostic. The new run uses the steerable
`relative_gluon` metric and the full x grid through x=1.

## High-x gluon closure

The central fitted HERAPDF was exported at `Q=1.65 GeV` through `x=1` and
projected without PDF uncertainty members. Reproduce the test with:

```bash
cd fits/SMP-22-014/figure2/likelihood_scans
python high_x_gluon_closure.py
```

Outputs and full numerical metrics are stored in `high_x_gluon_closure/`. The
comparison separates three questions:

1. The current 140-node `dist0` projection, ending at `x=0.83269`, reproduces
   its stored coefficients exactly but its gluon turns negative near `x=0.8`.
2. Extending the same absolute `dist0` metric through `x=1` improves but does
   not repair the relative high-x closure, because tiny absolute PDFs receive
   negligible weight.
3. The configurable `relative_gluon` metric adds a relative-gluon residual on `0.05<=x<=0.99`
   with weight 0.1. The same 100-mode basis then closes the gluon to at most
   0.82% on `0.6<=x<0.8`, 0.49% on `0.8<=x<0.9`, and 1.25% on
   `0.9<=x<0.99`, while changing the global residual/target-shift norm only
   from `4.53e-4` to `5.16e-4`.

This demonstrates that the basis can represent the fitted high-x gluon. The
relative-gluon term changes only how external PDFs are projected; it is not an
experimental likelihood contribution or a PDF prior.

The active configuration also has an optional relative light-valence term for
`u_v=u-ubar` and `d_v=d-dbar`. A three-point ADbar test is stored in
`adbar_gluon_valence_3point/`: it improves valence-PDF closure substantially,
but does not change the fixed-nuisance likelihood bias. This rules out alpha_s,
evolution amplification, and valence closure alone as the remaining dominant
source; the next metric candidate is the charge-weighted light-sea DIS
combination.

The next tested term is a relative photon-exchange `F2` proxy, the
charge-weighted sum `sum e_f^2 (q_f+qbar_f)` on `1e-4<x<0.1`. With weight 3 it
reduces the three-point ADbar HERA endpoint bias from about `+16.9/-18.5` to
`+13.4/-13.2`, and improves CMS at the same time. The plot is
`adbar_gluon_valence_f2_3point/plots/likelihood_ADbar.png`. This is useful but
not full closure; a future metric should use evolved, reaction-level HERA
structure functions rather than the input-scale proxy.
See `high_x_gluon_closure/high_x_gluon_closure.png` and
`high_x_gluon_closure/closure_metadata.yaml`.

## Files and interpretation

- `evaluate_direct.py`: immutable direct point record and exact-grid PDF export.
- `evaluate_pod.py`: immutable full-100-mode POD point record.
- `scan_tools.py`: likelihood parsing, grid parsing, checksums, and xFitter run helpers.
- `run_likelihood_scan.py`: projection, resumption, and compact aggregation.
- `plot_likelihood_scan.py`: one PNG per direction plus a multipage PDF.
- `high_x_gluon_closure.py`: central-PDF full-x gluon closure and metric diagnostic.
- `scan_results.npz`: numerical source for every plot; schema in `result_schema.yaml`.
- `scan_metadata.yaml`: fit ID, projection contract, software/config checksums, and status.

In the upper plot row both curves are shown relative to the direct likelihood at
the zero coordinate. The lower row is the pointwise bias diagnostic
`chi2_full_POD - chi2_direct`. A constant offset is therefore visible and is not
silently subtracted. HERA, CMS inclusive jets, and their sum are shown separately.

The reference point currently gives:

```text
direct:    1302.417647310198
full POD:  1301.164667751371
difference: -1.252979558827
```

This difference is a measured outcome of the projection/interface/evolution
test, not a closure assumption. The 100-mode projection's relative residual at
the reference is `3.2097e-4` in the configured metric.

The full POD solve is numerically sensitive (`cond(X) ~ 2.2e7`). For this
reason xFitter's LHAPDF6 writer was extended to accept explicit `Xvalues` and
`Qvalues` and to render IEEE-754 doubles with 17 digits. Do not replace this
with the ordinary `pdfs_q2val` text output, which keeps too little precision for
the high modes.
