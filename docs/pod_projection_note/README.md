# POD projection note

This is deliberately self-contained: equations, configuration choices, caveats,
and a compact bibliography are all in `main.tex`.

Compile from this directory:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Regenerate the likelihood-bias schematic with:

```bash
python plot_likelihood_bias_schematic.py
```
