#!/usr/bin/env bash
# Validate an already prepared checkout.  This script never downloads data,
# builds code, or runs a likelihood evaluation.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
with_pod=false
scan_input=""

usage() {
  cat <<'EOF'
Usage: ./scripts/verify-local.sh [--with-pod] [--plot SCAN_RESULTS.npz]

Checks the installed xFitter executable, the already downloaded public inputs,
Python syntax, and the deterministic Figure-2 projection-metric test.
--with-pod additionally validates an installed 250503_pod_basis_40k set.
--plot renders stored likelihood results into a temporary directory.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-pod) with_pod=true ;;
    --plot)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      scan_input="$2"
      shift
      ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
  shift
done

source "${project_root}/scripts/activate.sh"

for required in \
  "${project_root}/install/xfitter/bin/xfitter" \
  "${project_root}/xfitter-datafiles/hera/h1zeusCombined/inclusiveDis/1506.06042/HERA1+2_NCep_920-thexp.dat" \
  "${project_root}/xfitter-datafiles/lhc/cms/jets/2111.10431/Run2016_NNLO_y0.dat" \
  "${project_root}/xfitter-datafiles/lhc/cms/jets/2111.10431/FastNLO/1jet.NNLO.fnl5332h_y0_ptjet.tab"; do
  [[ -s "${required}" ]] || { echo "Missing required local input: ${required}" >&2; exit 1; }
done

head -c 64 "${project_root}/xfitter-datafiles/lhc/cms/jets/2111.10431/FastNLO/1jet.NNLO.fnl5332h_y0_ptjet.tab" \
  | grep -q 'version https://git-lfs.github.com/spec' \
  && { echo "The fastNLO grid is still a Git LFS pointer." >&2; exit 1; }

python -m py_compile \
  "${project_root}/pod_projection/"*.py \
  "${project_root}/fits/SMP-22-014/figure2/likelihood_scans/"*.py
(
  cd "${project_root}/fits/SMP-22-014/figure2/likelihood_scans"
  python test_projection_metrics.py
)

if "${with_pod}"; then
  "${project_root}/fits/SMP-22-014/figure2/pod_5mode/validate_basis.sh"
fi

if [[ -n "${scan_input}" ]]; then
  [[ -s "${scan_input}" ]] || { echo "No scan-result file: ${scan_input}" >&2; exit 1; }
  output_dir="$(mktemp -d "${TMPDIR:-/tmp}/incljets-likelihood-plot.XXXXXX")"
  python "${project_root}/fits/SMP-22-014/figure2/likelihood_scans/plot_likelihood_scan.py" \
    --input "${scan_input}" --output-dir "${output_dir}"
  echo "Plot validation output: ${output_dir}"
fi

echo "Local InclJets validation: PASS"
