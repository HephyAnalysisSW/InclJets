#!/usr/bin/env bash
# Source this file from a shell or an analysis wrapper after activating the
# Conda environment that provides ROOT, LHAPDF, GSL, yaml-cpp and a Fortran
# compiler.

_incljets_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export INCLJETS_ROOT="$(cd "${_incljets_script_dir}/.." && pwd)"
unset _incljets_script_dir

export PATH="${INCLJETS_ROOT}/install/xfitter/bin:${INCLJETS_ROOT}/install/qcdnum/bin:${PATH}"

_incljets_libs="${INCLJETS_ROOT}/install/xfitter/lib:${INCLJETS_ROOT}/install/qcdnum/lib"
if [[ -n "${CONDA_PREFIX:-}" ]]; then
  _incljets_libs="${CONDA_PREFIX}/lib:${_incljets_libs}"
fi

case "$(uname -s)" in
  Darwin) export DYLD_LIBRARY_PATH="${_incljets_libs}${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}" ;;
  *)      export LD_LIBRARY_PATH="${_incljets_libs}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" ;;
esac
unset _incljets_libs
