# Five-mode POD xFitter test

This is a prepared, deliberately unrun HERA+CMS fit using modes 1--5 of
`250503_pod_basis_40k`. The native POD definition and the external-projection
contract are in `pod.yaml`; fit/evolution choices are in `parameters.yaml`.

Safe PDF-level validation:

```bash
./validate_basis.sh
```

`run.sh` has an intentional guard. After the cards and coefficient steps have
been reviewed, the fit can be started explicitly with `./run.sh --run-fit`.

The validation first tests the native LHAPDF decomposition and then loads this
card through xFitter's normal plugin/YAML path. It runs QCDNUM evolution for
zero, one-hot, and mixed coefficient vectors at Q = 1.65, 10, 100, and 1000
GeV. It does not read datasets, construct a likelihood, or execute MINUIT
commands.

Change `basis_members` in `pod.yaml` to select modes. Omitting it, or setting
it to `all`, selects every non-reference member. `active_flavors: all` applies
the mode shifts to every flavor in this basis. `q0` is configurable but, by
default, must agree with the global xFitter `Q0` in `parameters.yaml`.

The POD member metadata uses `alpha_s(mZ)=0.118`, `mc=1.51 GeV`, and
`mb=4.92 GeV`. The analysis currently keeps the publication mass choices from
`../central/constants.yaml` (`mc=1.47 GeV`, `mb=4.50 GeV`) while fixing
`alpha_s(mZ)=0.118`. This distinction is intentional: the POD supplies the
input PDF at `Q0=1.65 GeV`, and QCDNUM applies the analysis evolution settings.
