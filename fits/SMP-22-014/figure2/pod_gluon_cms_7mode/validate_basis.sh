#!/usr/bin/env bash
set -euo pipefail

run_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$run_dir/../../../.." && pwd)"
source "$project_root/scripts/activate.sh"

"$project_root/install/xfitter/bin/pod-basis-check" \
  gluon_POD_nongluon_PDF4LHC21 1.65 0 1 2 3 4 5 6 7

validation_dir="$(mktemp -d "${TMPDIR:-/tmp}/xfitter-pod-validation.XXXXXX")"
cleanup() {
  case "$validation_dir" in
    */xfitter-pod-validation.*) rm -rf -- "$validation_dir" ;;
  esac
}
trap cleanup EXIT
cp -L "$run_dir/parameters.yaml" "$validation_dir/parameters.yaml"
cp -L "$run_dir/pod.yaml" "$validation_dir/pod.yaml"
cp -L "$run_dir/constants.yaml" "$validation_dir/constants.yaml"
cd "$validation_dir"
"$project_root/install/xfitter/bin/pod-xfitter-check"
