#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/project_pdf.py"
python "$script_dir/plot_ct18_pod_projection_uncertainties.py"
