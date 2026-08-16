# Inclusive jets with xFitter and POD PDFs

This repository reconstructs the NNLO HERA+CMS inclusive-jet PDF fit behind
Figure 2 of the CMS inclusive-jet addendum, JHEP 12 (2022) 035. It also adds a
POD PDF parameterisation to xFitter and tools for projecting external PDFs onto
that parameterisation.

## Repository layout

| Path | Contents |
|---|---|
| `xfitter/` | Pinned xFitter submodule, on the `pod-decomposition` fork branch. It contains the POD decomposition and supporting fixed-likelihood/output changes. |
| `pod_projection/` | Python projection utilities and public CT18 validation example. |
| `fits/SMP-22-014/figure2/` | CMS/HERA cards, frozen reference state, likelihood scripts, and documentation. |
| `xfitter-datafiles/` | Public xFitter datafiles checkout. It is downloaded locally, never vendored into this repository. |
| `build/`, `install/` | Local build products; never committed. |

## Setup

Create the supplied Conda environment, then clone the repository recursively
and run the public setup script:

```bash
conda env create -f environment.yml
conda activate incljets

git clone git@github.com:HephyAnalysisSW/InclJets.git
cd InclJets
./scripts/setup-public.sh
source ./scripts/activate.sh
```

The script initializes the pinned xFitter submodule at shallow depth; builds the bundled public
QCDNUM 18-00/00 source into `install/qcdnum`; configures and installs xFitter
into `install/xfitter`; and checks out `xfitter-datafiles` revision
`4ed3a5d46872df39c82ed10f3aa9356f382f3c41`. The data checkout contains the
public HERA/CMS cards, correlations, corrections, and theory grids. It uses a
sparse checkout of only the HERA I+II and CMS inclusive-jet subtrees required
here. Do not add it to Git. Git LFS is required for the NNLO grids; set
`JOBS=N` to change the default two build jobs.

`environment.yml` is the supported public environment. An existing compatible
environment (including the local `root` environment used during development)
also works, but it must provide the packages listed there.

The first lightweight validation after public setup is the Figure-2 smoke fit:

```bash
cd fits/SMP-22-014/figure2/smoke
./run.sh
./plot.sh
```

To validate a prepared checkout later, without downloading, rebuilding, or
running a likelihood evaluation, use:

```bash
./scripts/verify-local.sh
```

If the private analysis POD set is installed, add `--with-pod` to check its
member layout through the xFitter interface. To verify plotting from an
existing result without rerunning a scan, pass its compact results file:

```bash
./scripts/verify-local.sh --with-pod \
  --plot fits/SMP-22-014/figure2/likelihood_scans/production_16x21_gluon_valence_f2/scan_results.npz
```

## Running the Figure-2 fit

The standard reconstructed HERA+CMS NNLO fit is intentionally separate from
the smoke test and takes hours on a laptop. From the repository root, after
public setup, run:

```bash
source ./scripts/activate.sh
cd fits/SMP-22-014/figure2/central
./run.sh
```

This writes the central fit under `central/output/` and validates the MINUIT
and HESSE result. To produce its Hessian PDF members and plots after a
successful central fit, use the documented continuation in
[`errors/README.md`](fits/SMP-22-014/figure2/errors/README.md).

## Likelihood comparison plot

The direct-versus-POD comparison needs the installed
`250503_pod_basis_40k` LHAPDF set in addition to the public setup. It fixes all
PDF and nuisance parameters at each point: it does not minimize or profile.
Run the short three-point smoke comparison and make its plot with:

```bash
source ./scripts/activate.sh
cd fits/SMP-22-014/figure2/likelihood_scans
python run_likelihood_scan.py \
  --output smoke_Bg_3point_gluon_valence_f2 \
  --parameters Bg \
  --coordinates=-0.5,0,0.5
python plot_likelihood_scan.py \
  --input smoke_Bg_3point_gluon_valence_f2/scan_results.npz \
  --output-dir smoke_Bg_3point_gluon_valence_f2/plots
```

The resulting `likelihood_Bg.png` and `likelihood_comparison.pdf` compare the
direct and full-POD HERA, CMS, and summed chi-square terms. The full 16
parameter x 21 point scan is documented in the
[`likelihood-scan README`](fits/SMP-22-014/figure2/likelihood_scans/README.md).

### POD inputs

The current analysis basis is `250503_pod_basis_40k`: member 0 is the
full-flavour reference and members 1--100 are its POD directions at
Q0 = 1.65 GeV. Its construction is implemented by the public
[WMIN model](https://github.com/HEP-PBSP/wmin-model), and the public
[NNPOD runcards](https://github.com/comane/NNPOD-wiki) use this exact set name.
The verified installed version used here has 101
members and whole-directory SHA-256 fingerprint
`b3eeae40a8c753090b22beeaf3adb393b36beb54c56367b159426db30c17cfd2`.

Install with:

```bash
./scripts/install-pod-basis.sh private_pdfs/250503_pod_basis_40k.tar.gz
```

Verify an installation and the xFitter interface with:

```bash
cd fits/SMP-22-014/figure2/pod_5mode
./validate_basis.sh
```

`gluon_POD_nongluon_PDF4LHC21` is a separate, private gluon-only comparison
basis (member 0 plus 30 directions). It is intentionally not distributed or
required for the central full-flavour likelihood studies.

## Current status

- The central Figure 2 HERA I+II plus CMS inclusive-jet NNLO fit has been
  reconstructed and frozen under `fits/SMP-22-014/figure2/reference_fit/`.
- `xfitter/pdfdecomps/POD/` implements the POD input-scale decomposition.
- `pod_projection/` projects installed LHAPDF PDFs onto the native POD basis.
- `fits/SMP-22-014/figure2/likelihood_scans/` evaluates fixed-nuisance direct
  and full-POD likelihoods without minimization or nuisance profiling.
- `fits/SMP-22-014/figure2/cms_closure/` traces the full-POD compression
  defect through the CMS jet bins with every active CMS nuisance fixed to zero.

The generic POD package implements only basis loading, grids, and standard
projection metrics. The Figure-2 likelihood directory owns the optional
relative gluon, valence, and F2 closure objectives, so they cannot silently
affect a generic external-PDF projection.

With the currently configured relative-gluon+valence+F2 projection objective,
the fresh fixed-nuisance reference is total 1301.0235334369384 (HERA
1229.0831474048314, CMS 71.94038603210495). The earlier
relative-gluon-only reference, retained as a diagnostic, is
1301.1646677513709 (HERA 1229.7548014501342, CMS 71.40986630123827).

Detailed commands and validation notes live in:

- [`POD projection`](pod_projection/README.md)
- [`POD xFitter interface`](xfitter/pdfdecomps/POD/README.md)
- [`Figure 2 likelihood scans`](fits/SMP-22-014/figure2/likelihood_scans/README.md)
- [`CMS-only POD closure`](fits/SMP-22-014/figure2/cms_closure/README.md)

## Reproducibility policy

Source, configuration, compact frozen reference data, and checksums belong in
Git. Generated xFitter output directories, raw scans, derived plots/build
trees, public-data clones, and private PDF sets do not. `scripts/setup-public.sh`
and the smoke fit provide the portable public setup and validation path.

## References

- CMS measurement/addendum: [JHEP 12 (2022) 035](https://doi.org/10.1007/JHEP12(2022)035)
- CMS public result: [SMP-20-011](https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-20-011/)
- HERA I+II inclusive DIS combination: [arXiv:1506.06042](https://arxiv.org/abs/1506.06042)
- xFitter: <https://gitlab.cern.ch/fitters/xfitter>
