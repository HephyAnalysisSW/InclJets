#!/usr/bin/env bash
set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

test -s output/minuit.out.txt
test -s output/Results.txt
test -s output/fittedresults.txt
test -s run.log

grep -Eq 'FROM HESSE +STATUS=OK' output/minuit.out.txt
grep -Eq 'Covariance matrix status = +3 +16' run.log

member_count=$(find output -maxdepth 1 -type f -name 'pdfs_q2val_s??s_*.txt' -size +0 | wc -l | tr -d ' ')
if [ "${member_count}" -ne 80 ]; then
  echo "Expected 80 nonempty shifted PDF tables (16 members x 5 Q2 values); found ${member_count}." >&2
  exit 1
fi

echo "OK: one accurate 16-parameter HESSE and 16 symmetric PDF error members."
