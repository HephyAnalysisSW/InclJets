# POD projection note

This is deliberately self-contained: equations, configuration choices, caveats,
a code-oriented workflow map, and a compact bibliography are all in `main.tex`.
The PDF contains clickable source links. The corresponding repository paths are
also available directly:

- [Central analytic PDF-fit parameters](../../fits/SMP-22-014/figure2/central/parameters.yaml), [constants](../../fits/SMP-22-014/figure2/central/constants.yaml), and [steering card](../../fits/SMP-22-014/figure2/central/steering.txt)
- [Frozen reference fit and nuisances](../../fits/SMP-22-014/figure2/reference_fit/)
- [Projection configuration](../../fits/SMP-22-014/figure2/likelihood_scans/scan_config.yaml), [CMS metric](../../fits/SMP-22-014/figure2/likelihood_scans/projection_metrics.py), and [scan drivers](../../fits/SMP-22-014/figure2/likelihood_scans/)
- [Generic POD projection](../../pod_projection/pod_projection.py), [POD xFitter implementation](../../xfitter/pdfdecomps/POD/), and [POD-mediated analytic fit driver](../../fits/SMP-22-014/figure2/pod_mediated_fit/run_fit.py)

Compile from this directory:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Regenerate the likelihood-bias schematic with:

```bash
python plot_likelihood_bias_schematic.py
```
