#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00

# Post-MD structural clustering (CLoNe + PCA). Submitted with an afterok
# dependency on all production repeats, so it runs once MD and the post-processing
# that builds the combined/aligned trajectory have completed.
#
# It delegates to the Python package, sweeping the CLoNe density parameter
# (pdc = 1, 3, 5, 7) on the C-alpha coordinates after a 2-component PCA, matching
# the legacy run_structural.py invocation. Each pdc value writes its own results
# directory.

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

# The combined trajectory + topology produced by post-processing (trjcat/align).
# Overridable via CLUSTER_TRAJ / CLUSTER_TOPO in configurations.txt.
TRAJ="${CLUSTER_TRAJ:-$SIMULATION_DIR/production/combined_aligned.xtc}"
TOPO="${CLUSTER_TOPO:-$SIMULATION_DIR/production/final.pdb}"

if [ ! -f "$TRAJ" ] || [ ! -f "$TOPO" ]; then
    echo "Combined trajectory ($TRAJ) or topology ($TOPO) not found; skipping clustering."
    exit 0
fi

for pdc in 1 3 5 7; do
    out_dir="$SIMULATION_DIR/production/clustering/pdc_${pdc}"
    PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
      -m automd_saxs cluster \
      --traj "$TRAJ" \
      --topo "$TOPO" \
      --out-dir "$out_dir" \
      --at-sel "name CA" \
      --pca 2 \
      --pdc "$pdc" \
      || echo "clustering with pdc=${pdc} failed; continuing"
done
