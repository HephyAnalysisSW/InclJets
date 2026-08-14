#!/usr/bin/env bash
# Install only public source/data under this checkout.  Private POD LHAPDF
# grids are intentionally out of scope; see the top-level README.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
jobs="${JOBS:-2}"

for tool in cmake git make tar; do
  command -v "${tool}" >/dev/null || {
    echo "Missing required command: ${tool}" >&2
    exit 1
  }
done
git lfs version >/dev/null || {
  echo "Git LFS is required to download the public NNLO fastNLO grids." >&2
  exit 1
}

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "Activate the Conda environment with ROOT/LHAPDF/GSL/yaml-cpp first." >&2
  exit 1
fi

git -C "${project_root}" submodule update --init --recursive

data_revision="4ed3a5d46872df39c82ed10f3aa9356f382f3c41"
data_dir="${project_root}/xfitter-datafiles"
if [[ ! -d "${data_dir}/.git" ]]; then
  git clone --filter=blob:none https://gitlab.cern.ch/fitters/xfitter-datafiles.git "${data_dir}"
fi
git -C "${data_dir}" fetch --depth=1 origin "${data_revision}"
git -C "${data_dir}" checkout --detach "${data_revision}"
git -C "${data_dir}" lfs pull

qcdnum_prefix="${project_root}/install/qcdnum"
if [[ ! -x "${qcdnum_prefix}/bin/qcdnum-config" ]]; then
  qcdnum_build="${project_root}/build/qcdnum"
  mkdir -p "${qcdnum_build}"
  qcdnum_source="${qcdnum_build}/qcdnum-18-00-00"
  if [[ ! -d "${qcdnum_source}" ]]; then
    tar -xzf "${project_root}/xfitter/tools/qcdnum180000.tar.gz" -C "${qcdnum_build}"
  fi
  (
    cd "${qcdnum_source}"
    ./configure --prefix="${qcdnum_prefix}"
    make -j "${jobs}" install
  )
fi

export PATH="${qcdnum_prefix}/bin:${PATH}"
cmake -S "${project_root}/xfitter" -B "${project_root}/build/xfitter" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${project_root}/install/xfitter" \
  -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
  -Dyaml-cpp_DIR="${CONDA_PREFIX}"
cmake --build "${project_root}/build/xfitter" --parallel "${jobs}"
cmake --install "${project_root}/build/xfitter"

echo "Public setup complete. In a new shell: source ${project_root}/scripts/activate.sh"
