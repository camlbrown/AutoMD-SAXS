#!/bin/bash

#SBATCH --ntasks-per-node=12
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=scavenger


SIMULATION_DIR="$1"
config_file="$SIMULATION_DIR/configurations.txt"
if [ -f "$config_file" ]; then
    source "$config_file"
else
    echo "Error: configurations.txt file not found in $SIMULATION_DIR."
    exit 1
fi

GMXLIB="$BASE_DIR/ff_files"
export GMXLIB

load_gmx

cd $PRODUCTION_DIR

cp $NPT_DIR/npt.gro .  && cp $NPT_DIR/topol.top . && cp $NPT_DIR/npt.cpt . && cp $NPT_DIR/*itp . && cp $NPT_DIR/*txt . && cp $NPT_DIR/*ndx .

gmx_mpi grompp -f $MDP_DIR/lig_md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 100 2>&1 | tee md_grompp.log

cp md.tpr $REPEAT_DIR1 && cp md.tpr $REPEAT_DIR2 && cp md.tpr $REPEAT_DIR3

cp $NPT_DIR/topol.top $REPEAT_DIR1 && cp $NPT_DIR/topol.top $REPEAT_DIR2 && cp $NPT_DIR/topol.top $REPEAT_DIR3 
