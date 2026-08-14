# Inclusive jets with xFitter and POD PDFs

This repository reconstructs the NNLO HERA+CMS inclusive-jet PDF fit behind
Figure 2 of the CMS inclusive-jet addendum, JHEP 12 (2022) 035. It also adds a
POD PDF parameterisation to xFitter and tools for projecting external PDFs onto
that parameterisation.

The repository is being prepared for sharing. The setup below covers the
public inputs; the private POD LHAPDF sets are intentionally not distributed.

## Repository layout

| Path | Contents |
|---|---|
| `xfitter/` | xFitter source with the POD decomposition and small supporting changes. This will be maintained as a pinned fork/submodule. |
| `pod_projection/` | Python projection utilities and public CT18 validation example. |
| `fits/SMP-22-014/figure2/` | CMS/HERA cards, frozen reference state, likelihood scripts, and documentation. |
| `xfitter-datafiles/` | Public xFitter datafiles checkout. It is downloaded locally, never vendored into this repository. |
| `build/`, `install/` | Local build products; never committed. |

## Setup

Clone this repository, then fetch the public xFitter input data at the pinned
revision used for the analysis:

```bash
git clone git@github.com:HephyAnalysisSW/InclJets.git
cd InclJets

git clone --filter=blob:none \
  https://gitlab.cern.ch/fitters/xfitter-datafiles.git xfitter-datafiles
git -C xfitter-datafiles checkout 4ed3a5d46872df39c82ed10f3aa9356f382f3c41
```

The data checkout contains the public HERA and CMS cards, correlations,
corrections, and theory grids. Do not add it to Git.

### Private POD inputs

The following LHAPDF sets are required for POD studies but are private and are
not included here:

- `250503_pod_basis_40k` — full-flavour reference plus 100 POD directions;
- `gluon_POD_nongluon_PDF4LHC21` — gluon-only reference plus 30 directions.

Install them through an authorized local LHAPDF path. The code checks the set
name and expected member range at runtime. A checksummed distribution procedure
will be documented here once redistribution/access is agreed.

## Current status

- The central Figure 2 HERA I+II plus CMS inclusive-jet NNLO fit has been
  reconstructed and frozen under `fits/SMP-22-014/figure2/reference_fit/`.
- `xfitter/pdfdecomps/POD/` implements the POD input-scale decomposition.
- `pod_projection/` projects installed LHAPDF PDFs onto the native POD basis.
- `fits/SMP-22-014/figure2/likelihood_scans/` evaluates fixed-nuisance direct
  and full-POD likelihoods without minimization or nuisance profiling.

Detailed commands and validation notes live in:

- [`POD projection`](pod_projection/README.md)
- [`POD xFitter interface`](xfitter/pdfdecomps/POD/README.md)
- [`Figure 2 likelihood scans`](fits/SMP-22-014/figure2/likelihood_scans/README.md)

## Reproducibility policy

Source, configuration, compact frozen reference data, checksums, and selected
plots belong in Git. Generated xFitter output directories, raw scans, build
trees, public-data clones, and private PDF sets do not. A small smoke test and
a portable build recipe are the next repository milestones.

## References

- CMS measurement/addendum: [JHEP 12 (2022) 035](https://doi.org/10.1007/JHEP12(2022)035)
- CMS public result: [SMP-20-011](https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-20-011/)
- HERA I+II inclusive DIS combination: [arXiv:1506.06042](https://arxiv.org/abs/1506.06042)
- xFitter: <https://gitlab.cern.ch/fitters/xfitter>
