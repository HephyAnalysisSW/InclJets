# Figure 2 Hessian error-member continuation

This directory reproducibly generates Hessian PDF error members from the converged central minimum without repeating the multi-hour minimisation.

## Why this is separate

The first central card contained both an explicit `hesse` command and `doErrors: Hesse`. In xFitter 2.2.1, `doErrors: Hesse` invokes HESSE internally before diagonalising the covariance and producing shifted PDFs. The original job therefore ran HESSE twice after MIGRAD. The first pass was accurate; the redundant second pass inherited much smaller finite-difference steps, became non-positive-definite, and stopped before writing error members.

This continuation embeds the converged values from `../central/output/parsout_1`, retains the original robust MINUIT step sizes, and deliberately contains no explicit `hesse` command. `call fcn 3` stores the central prediction, after which `doErrors: Hesse` performs the sole HESSE and constructs 16 symmetric eigenvector members.

No data, theory, chi2, PDF-form, mass, scale, correction, or evolution setting differs from the central fit.

## Re-run

```sh
cd fits/SMP-22-014/figure2/errors
./run.sh
```

`run.sh` records the full console stream in `run.log` and then invokes
`check.sh`. The check requires an accurate 16-parameter covariance and 80
nonempty shifted tables: 16 members at each of the five requested Q2 values.
`./check.sh` can also be rerun independently without repeating the fit.

Generate post-fit PDF plots with the symmetric Hessian uncertainty band using:

```sh
./plot.sh
```

The plotting wrapper reads the central and 16 shifted PDF tables in `output/`
and writes a combined `plots/plots.pdf` plus individual PNG files for every
stored Q2 value. It does not rerun xFitter or the Hessian calculation.

xFitter writes `Status.out` before its separate error-analysis phase. Because this continuation intentionally does not run MIGRAD, that file may say `Failed` even when the subsequent HESSE and member generation succeed. The authoritative continuation checks are the HESSE status in `output/minuit.out.txt`, `Covariance matrix status = 3 16` in `run.log`, and the shifted-file count enforced by `check.sh`.

For future full fits, use only `migrad ...` and `call fcn 3` in `MINUIT.Commands`; leave `doErrors: Hesse` enabled and do not also add an explicit `hesse` command.

## Validated result (2026-08-05)

The continuation finished with exit code 0. The sole post-fit HESSE used 245
function calls (246 including the initial central FCN), reported
`STATUS=OK`, EDM = 4.1e-5, and covariance status 3 for all 16 parameters.
It then wrote 16 symmetric shifted members, each at five Q2 values. The final
displayed evaluation was 278 and the measured runtime was about 684 seconds
(11.4 minutes). The recovered result is

```text
chi2/Ndof = 1302.418/1118
alpha_s(mZ) = 0.116638 +/- 0.001549
```
