#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=00:30:00

# Post-MD SAXS analysis. Submitted by run_MD.sh / `automd_saxs run` with an
# afterok dependency on all production repeats, so it runs automatically once MD
# (and the post-processing that produces the crysol fit outputs) has completed.
#
# It delegates to the Python package, which collects crysol_summary.txt / *.fit /
# Rg_distr*.txt from the SAXS directory, writes a chi^2-ranked saxs_summary.csv,
# and folds bestChi2 / bestFrame / rgMean into an analysis manifest.

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <simulation_directory>"
    exit 1
fi

SIMULATION_DIR="$1"
config_file="$SIMULATION_DIR/configurations.txt"
if [ ! -f "$config_file" ]; then
    echo "Error: configurations.txt not found in $SIMULATION_DIR" >&2
    exit 1
fi
source "$config_file"

if [ "${SAXS_FILE:-None}" = "None" ]; then
    echo "No SAXS data for this job; skipping analysis."
    exit 0
fi

shopt -s nullglob
dats=("$SAXS_DIR"/*.dat)
shopt -u nullglob
SAXS_ARG=()
[ "${#dats[@]}" -gt 0 ] && SAXS_ARG=(--saxs "${dats[0]}")

PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
  -m automd_saxs analyze \
  --analysis-dir "$SAXS_DIR" \
  --out-dir "$SAXS_DIR" \
  --config "$config_file" \
  "${SAXS_ARG[@]}"
