#!/usr/bin/env bash
# Install a privately shared LHAPDF POD archive without placing it in Git.
set -euo pipefail

set_name="250503_pod_basis_40k"
expected_sha256="c26b193ae8355fd8e3a305f0b806704b0583609324685655553b3c0689f1b21d"
archive=""
destination=""

usage() {
  cat <<EOF
Usage: ./scripts/install-pod-basis.sh ARCHIVE [--destination LHAPDF_PATH] [--sha256 SHA256]

Installs the privately shared ${set_name} LHAPDF archive. By default the
archive must have this SHA-256:
  ${expected_sha256}

The default destination is the first LHAPDF search path of the active Python
environment. The installer refuses to overwrite an existing set.
EOF
}

[[ $# -ge 1 ]] || { usage >&2; exit 2; }
archive="$1"
shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --destination)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      destination="$2"
      shift
      ;;
    --sha256)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      expected_sha256="$2"
      shift
      ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
  shift
done

[[ -f "$archive" ]] || { echo "Archive not found: $archive" >&2; exit 1; }
if [[ -z "$destination" ]]; then
  destination="$(python - <<'PY'
import lhapdf
print(lhapdf.paths()[0])
PY
)"
fi
[[ -d "$destination" ]] || { echo "LHAPDF destination does not exist: $destination" >&2; exit 1; }
[[ ! -e "$destination/$set_name" ]] || { echo "Refusing to overwrite: $destination/$set_name" >&2; exit 1; }

if command -v shasum >/dev/null; then
  actual_sha256="$(shasum -a 256 "$archive" | awk '{print $1}')"
else
  actual_sha256="$(sha256sum "$archive" | awk '{print $1}')"
fi
[[ "$actual_sha256" == "$expected_sha256" ]] || {
  echo "Archive SHA-256 mismatch: got $actual_sha256" >&2
  exit 1
}

tar -tzf "$archive" | awk -v root="$set_name/" '
  $0 !~ /^250503_pod_basis_40k\// || $0 ~ /(^|\/)\.\.($|\/)/ { exit 1 }
  END { if (NR == 0) exit 1 }
' || { echo "Archive has unsafe or unexpected paths." >&2; exit 1; }

tar -xzf "$archive" -C "$destination"
member_count="$(find "$destination/$set_name" -maxdepth 1 -type f -name "${set_name}_*.dat" | wc -l | tr -d ' ')"
[[ "$member_count" == 101 ]] || {
  echo "Expected 101 PDF members after extraction; found $member_count." >&2
  exit 1
}
[[ -s "$destination/$set_name/${set_name}.info" ]] || {
  echo "Missing LHAPDF metadata file after extraction." >&2
  exit 1
}
echo "Installed $set_name (101 members) in $destination/$set_name"
