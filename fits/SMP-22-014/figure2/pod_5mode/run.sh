#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" != "--run-fit" ]]; then
  echo "Fit intentionally not started. Review the POD cards, then run:"
  echo "  ./run.sh --run-fit"
  exit 2
fi

run_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${run_dir}/../../../.." && pwd)"
source "${project_root}/scripts/activate.sh"
cd "$run_dir"

xfitter 2>&1 | tee run.log
