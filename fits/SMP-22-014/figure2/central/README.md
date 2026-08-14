# SMP-22-014 Figure 2: central NNLO fit

This directory is a reconstruction of the central HERA+CMS fit used for
Figure 2 of the addendum to CMS-SMP-20-011 (JHEP 12 (2022) 035). It is kept
outside the xFitter source checkout so the analysis configuration can be
versioned independently of the fitter.

The public repositories provide all experimental inputs, statistical
correlation matrices, NP/EW corrections, and the NNLO fastNLO grids. They do
not provide the publication's master `steering.txt` and `parameters.yaml`.
Consequently, the following choices in these cards are taken directly from the
paper:

- HERA I+II inclusive DIS plus the four CMS R=0.7 inclusive-jet rapidity bins;
- NNLO evolution and NNLO fastNLO interpolation grids, with jet pT as the
  central renormalisation and factorisation scale;
- Q0^2 = 1.9 GeV^2, Qmin^2 = 7.5 GeV^2, mc = 1.47 GeV, mb = 4.5 GeV,
  and fs = 0.4;
- the PDF forms in Eqs. (11)-(15), giving 15 free PDF-shape parameters, plus
  fitted alpha_s(mZ);
- Hessian fit uncertainties with Delta chi^2 = 1.

The numerical starting values and MINUIT command sequence are reconstructed
from current xFitter HERA/jet examples. They are seeds, not published fit
results. The chi-square settings are the standard multiplicative/Hessian
settings used by xFitter's public HERA examples; these should be confirmed
against the original production cards if those become available.

## Current state

The `datafiles` link resolves to the sparse checkout at the project root. All
text inputs and the four SHA-256-verified NNLO grids are present. Native ARM64
builds are installed at `../../../../install/xfitter` and
`../../../../install/qcdnum`.

The first central fit completed on 2026-08-05 with chi2/Ndof = 1302.42/1118
and alpha_s(mZ) = 0.116638. MIGRAD converged after 2651 total function calls,
and its first post-fit HESSE was accurate. The original card nevertheless
contained both an explicit `hesse` command and `doErrors: Hesse`; the latter
runs HESSE internally, so a redundant second pass became non-positive-definite
and did not write shifted PDF members. The central values in `output/` remain
valid. The successful continuation in `../errors` recovered an accurate
16-parameter covariance, alpha_s(mZ) = 0.116638 +/- 0.001549, and all 16
symmetric PDF error members.

## Reproducible run sequence

For a new full central fit, run:

```sh
cd fits/SMP-22-014/figure2/central
./run.sh
```

The active `parameters.yaml` contains the corrected MINUIT procedure:

```text
set str 2
migrad 200000
call fcn 3
doErrors: Hesse
```

Here `doErrors: Hesse` is a xFitter directive, not a MINUIT command: after the
command block completes it invokes the one required post-fit HESSE pass and
constructs the shifted PDF members. Do not add another explicit `hesse` to the
command block. `run.sh` records the console stream in `run.log` and calls
`check.sh`, which requires converged MIGRAD, a successful final HESSE, an
accurate 16-parameter covariance, and 80 nonempty shifted tables (16 members
at five requested Q2 values).

To regenerate only the error members from the already converged minimum, use
the much shorter `../errors/run.sh`; its seeds are frozen to the central
`parsout_1` values and its procedure is described in `../errors/README.md`.
The 2026-08-05 continuation took about 11.4 minutes and ended at displayed
evaluation 278, compared with several hours for the full minimisation.

After either workflow has produced validated error members, make the post-fit
PDF plots from the continuation directory:

```sh
cd fits/SMP-22-014/figure2/errors
./plot.sh
```

This reads the central and 16 symmetric Hessian-member tables without rerunning
the fit. It writes the multipage `plots/plots.pdf` and individual PNG plots at
all five stored Q2 values. The wrapper options are kept in `plot.sh` so plot
formatting can be changed and regenerated independently of the fit.

## Physics and software chain

At the starting scale Q0^2 = 1.9 GeV^2, xFitter constructs the gluon, valence,
and sea distributions from the forms in Eqs. (11)-(15) of the CMS paper.
Valence-number and momentum sum rules determine the dependent normalisations;
15 PDF-shape parameters and alpha_s(mZ) are varied. For every MINUIT function
call, QCDNUM evolves the trial PDFs and alpha_s at NNLO and xFitter evaluates
all DIS and jet predictions, constructs the correlated chi-square, and returns
its scalar value to MINUIT. This repeated full theory evaluation explains the
runtime of the fit.

The seven HERA cards are the combined H1+ZEUS inclusive DIS measurement from
arXiv:1506.06042: NC e+p at proton energies 920, 820, 575, and 460 GeV, NC e-p,
CC e+p, and CC e-p. Their cards contain the combined statistical, uncorrelated,
correlated, and procedural uncertainties. The Q2 > 7.5 GeV^2 cut leaves 1056
HERA points. Neutral-current reduced cross sections are routed through
`RT_DISNC` (the RT general-mass variable-flavour-number scheme selected by the
card), while charged-current cross sections use `BaseDISCC`; both consume the
NNLO PDFs evolved by QCDNUM.

The CMS input is the 2016 double-differential inclusive-jet measurement from
arXiv:2111.10431, using anti-kT jets with R = 0.7 and 33.5 fb^-1. The four data
cards cover |y| = 0.0-0.5, 0.5-1.0, 1.0-1.5, and 1.5-2.0 and contribute
22, 21, 19, and 16 active points. Their columns hold the measured cross
sections and experimental sources such as luminosity, JER, pileup, unfolding,
and the decomposed JES uncertainties. The ten `Run2016_yi_yj.corr` files are
the unique diagonal and lower-triangular blocks of the symmetric statistical
correlation matrix, including correlations between rapidity bins.

For each rapidity bin, the data card defines the theory expression

```
F * NP * EW
```

`F` is evaluated by the fastNLO 2.5.0 toolkit from one of these NNLO
interpolation tables:

- `FastNLO/1jet.NNLO.fnl5332h_y0_ptjet.tab`
- `FastNLO/1jet.NNLO.fnl5332h_y1_ptjet.tab`
- `FastNLO/1jet.NNLO.fnl5332h_y2_ptjet.tab`
- `FastNLO/1jet.NNLO.fnl5332h_y3_ptjet.tab`

The tables encode the NNLO partonic calculation and fastNLO convolves them
with the current trial PDFs and alpha_s using muR = muF = individual-jet pT.
The separate `NP/NP_y*.dat` factors correct the fixed-order prediction for
hadronisation and multiparton interactions; `EW/EW_y*.dat` applies the NLO
electroweak factor. The public checkout also contains NLO+NLL, contact-
interaction, and scale-varied cards, but none of them is used by this central
NNLO run.

The configured chi-square uses Poisson statistical scaling, linear
uncorrelated and correlated systematic scaling, the full CMS statistical
covariance, and Hessian nuisance parameters for correlated sources. MINUIT
96.03 runs strategy 2 on the 16 free parameters and stores the converged
prediction with a final FCN call. xFitter's `doErrors: Hesse` phase then runs
HESSE and constructs the symmetric shifted members. OpenBLAS/LAPACK provide
the covariance and Hessian linear algebra.
ROOT is linked and `xfitter-draw` is installed for later result inspection and
plots; it is not the source of the NNLO cross-section prediction.

The software versions in this local reproduction are xFitter 2.2.1, bundled
fastNLO 2.5.0_2826, and QCDNUM 18-00/00. The paper reports xFitter 2.2.1 with
QCDNUM 17-01/14. QCDNUM therefore remains a controlled implementation
difference to assess when comparing the converged chi-square and fitted
parameters.

All experimental and theory-table inputs were taken from the public
`xfitter-datafiles` repository. The jet checkout includes the repository update
labelled "Updating to NNLO grids to correspond to 2111.10431 addendum"
(`e7933d99...`); the four large table objects were independently SHA-256
verified after download.

## Changing the analysis

- Edit `steering.txt` to add/remove datasets or correlation matrices, change
  the chi-square prescription or cuts, and select the output Q2 and x grid.
- Edit `parameters.yaml` to free/fix PDF or alpha_s parameters, change their
  seeds and step sizes, alter PDF functional forms, select evolution/order, or
  change the MINUIT command sequence. When `doErrors: Hesse` is enabled, do
  not also put an explicit `hesse` in `MINUIT.Commands`.
- Edit `constants.yaml` for heavy-quark masses and electroweak/CKM constants.
- Public data cards, correlation matrices, NP/EW corrections, and fastNLO
  grids live below the `datafiles` symlink; avoid editing that checkout for a
  fit variation unless the input itself is intentionally being corrected.
- Copy this `central` directory to a sibling directory for each model, scale,
  parametrisation, or data variation. Keep its cards and `output/` together so
  results can always be traced back to the exact configuration.

The Apple-Silicon build uses compile-time capacities of 3000 data points, 1000
systematics, 50 datasets, and 20 chi2-scan points. Increase the corresponding
`NTOT_C`, `NSYSMAX_C`, `NSET_C`, or `NCHI2POINTS_C` CMake compiler definitions
and rebuild if a future analysis exceeds one of those limits.

## First validation targets

After the dependencies and grids are available, the central run should be
checked against Table B.1 of the addendum:

- global chi2 / Ndof: 1302 / 1118;
- CMS rapidity bins: 8.6 / 22, 23 / 21, 13 / 19, and 14 / 16;
- fitted alpha_s(mZ): 0.1166 with fit uncertainty 0.0014.

Agreement at this stage validates the central fit and Hessian band. The model,
scale, and parametrisation bands in Figure 2 require separate variation fits
and envelope construction.

## Sources

- Combined paper and addendum: https://arxiv.org/abs/2111.10431
- Addendum DOI: https://doi.org/10.1007/JHEP12(2022)035
- CMS public result: https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-20-011/
- HEPData record: https://www.hepdata.net/record/ins1972986
