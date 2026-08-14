#!/usr/bin/env bash
set -euo pipefail

run_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${run_dir}/../../../.." && pwd)"
source "${project_root}/scripts/activate.sh"
cd "${run_dir}"

exec xfitter-draw \
  --outdir plots \
  --bands \
  --filledbands \
  --q2all \
  --splitplots-png \
  --highres \
  --no-data \
  --no-shifts \
  --no-tables \
  output:'CMS inclusive jets NNLO fit'
