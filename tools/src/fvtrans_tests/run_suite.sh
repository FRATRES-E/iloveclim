#!/usr/bin/env bash
#-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#  Runs the FV transport evaluation suite (fvt_suite) and the plots (ROADMAP step 1.1, wiso:1c:0001).
#
#  Usage: ./run_suite.sh <outdir> [--full]
#    T21 (32 x 64) with 72 steps per period (the ECBilt 4 h step), winds at mid-step and at the start of the step,
#    T21 with 144 steps (half step), T42 with 144 and T85 with 288 steps (same Courant numbers as T21/72);
#    --full adds T170 with 576 steps (about one hour; T85 takes about 4 min, the rest seconds).
#  Outputs in <outdir>: summary.csv, one netCDF file per case/mode/scheme/run, figures in <outdir>/figures.
#-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "usage: $0 <outdir> [--full]"
  exit 2
fi
outdir=$1
here=$(cd "$(dirname "$0")" && pwd)

mkdir -p "${outdir}"
rm -f "${outdir}/summary.csv"
make -C "${here}" suite > /dev/null

runs=("32 64 72 mid" "32 64 72 start" "32 64 144 mid" "64 128 144 mid" "128 256 288 mid")
if [ "${2:-}" = "--full" ]; then
  runs+=("256 512 576 mid")
fi

for r in "${runs[@]}"; do
  set -- ${r}
  echo "fvt_suite $*"
  "${here}/fvt_suite" "$1" "$2" "$3" "$4" "${outdir}" > "${outdir}/log_n$1_s$3_$4.txt"
  tail -1 "${outdir}/log_n$1_s$3_$4.txt"
done

python3 "${here}/plot_suite.py" "${outdir}"
