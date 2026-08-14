# Immutable global-fit reference

Fit ID: `smp22014_global_hera_cms_20260805`

This directory freezes the completed HERA I+II plus CMS inclusive-jet global
minimum and its successful 16-parameter Hesse continuation. It is the sole
reference point for the no-minimization likelihood scans.

- `fit_result.yaml`: canonical parameters, fit summary, provenance, and hashes.
- `covariance.npz`: ordered values, Hesse errors, covariance, correlation, and eigenvalues.
- `nuisances.yaml`: all 198 global-best-fit experimental nuisance values as
  retained by xFitter's four-decimal `Results.txt` output.
- `cards/`: exact minimization and covariance-continuation cards.
- `raw/`: minimal xFitter outputs from which the snapshot was extracted.
- `files.sha256`: integrity hashes for every snapshot payload file.

Do not edit or regenerate this directory in place. The generator refuses to
overwrite it; create a new output directory and fit ID for a new reference.
No likelihood evaluation or minimization is performed by the generator.
The exact runtime nuisance array was not saved by the completed fit. Before a
scan uses the rounded values as fixed inputs, its evaluator must reproduce the
direct-reference chi2 within an explicitly recorded tolerance.
