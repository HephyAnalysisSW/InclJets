# POD-mediated analytic fit

The 16 HERAPDF analytic parameters are the only outer coordinates. At every
trial point, xFitter exports the analytic PDF at 1.65 GeV, Python projects it
onto the complete POD basis, and a second xFitter evaluation profiles the
HERA+CMS experimental nuisances. POD coefficients are derived, not fitted.

Prepare the central point:

    python run_fit.py --output nominal_full_pod --prepare-only

Start the outer fit:

    python run_fit.py --output nominal_full_pod --maxiter 200

Completed evaluations are recorded in history.json and a restarted command
begins from the best retained point. One point costs about 45 seconds locally.
