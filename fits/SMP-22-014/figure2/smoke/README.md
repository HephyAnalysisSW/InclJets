# Figure 2 smoke fit

This is a small, end-to-end debugging configuration derived from the reconstructed central fit. It is intended for testing xFitter changes and developing output/plotting scripts, not for physics comparisons with the publication.

It retains the important software path:

- 69 HERA NC e+p points with Q2 > 1000 GeV2, calculated with QCDNUM and RT_DISNC;
- all 22 CMS points in the first rapidity bin, including its 22-by-22 statistical correlation matrix;
- the y0 NNLO fastNLO grid and the same NP and electroweak correction tables as the central fit;
- the same chi2 prescription and a MINUIT/HESSE fit, but varying only Bg and alpha_s.

The directory reuses the central QCDNUM weight cache through `unpolarised.wgt` and links the common public inputs through `datafiles`. Its own `output/` is independent.

Run it with:

```sh
./run.sh
```

The validated run has 91 fitted points and takes about 82 seconds on this machine (68 FCN evaluations). The reduced parameterisation and data selection mean its fitted values and uncertainties have no publication-level interpretation.

Plot the fit with:

```sh
./plot.sh
```

This writes overview and individual PNG plots under `plots/`. For custom plot scripts, the most useful fixtures are `output/fittedresults.txt` for data/theory and `output/pdfs_q2val*.txt` for the central and two Hessian eigenvector PDF sets. `output/Results.txt` and `output/minuit.out.txt` contain the chi2 breakdown and fitted parameters.

Useful modifications:

- change the HERA threshold in `steering.txt` to trade speed for coverage;
- add another CMS `Run2016_NNLO_y*.dat` card and its correlation blocks to exercise multi-rapidity plotting;
- set a parameter's second array element in `parameters.yaml` to a nonzero MINUIT step to float it;
- omit `doErrors: Hesse` for the fastest central-value-only debugging cycle.
