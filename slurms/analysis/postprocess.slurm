#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --time=04:00:00

# Post-MD trajectory processing. Submitted with an afterok dependency on all
# production repeats. For each repeat it strips solvent, removes PBC jumps, fits
# out rotation/translation, builds a reference final.pdb, and extracts frames;
# then it concatenates the per-repeat final.xtc and aligns them into
# combined_aligned.xtc -- the input the cluster (and SAXS) stages consume.
#
# Delegates to the Python package, which builds the GROMACS command plan from the
# job layout and runs it. load_gmx makes gmx_mpi available on the compute node.

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
export GMXLIB="$BASE_DIR/ff_files"

# Make gmx_mpi available (no-op when GROMACS is a local install).
if command -v load_gmx >/dev/null 2>&1; then load_gmx; fi

PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
  -m automd_saxs postprocess \
  --config "$config_file" \
  --work-dir "$BASE_DIR"

if command -v unload_gmx >/dev/null 2>&1; then unload_gmx; fi
