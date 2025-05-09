#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=01:59:00
#SBATCH --partition=scavenger


# Source the configurations.txt file in the specified simulation directory
SIMULATION_DIR="$1"
config_file="$SIMULATION_DIR/configurations.txt"
if [ -f "$config_file" ]; then
    source "$config_file"
else
    echo "Error: configurations.txt file not found in $SIMULATION_DIR."
    exit 1
fi

load_gmx

cd $MINIM1_DIR && cp $GENION_DIR/*.gro . && cp $GENION_DIR/*.top .

gmx_mpi grompp -f $MDP_DIR/minim.mdp -c *.gro -p *.top -o minim1.tpr -maxwarn 5 2>&1 | tee minim1_grompp.log

gmx_mpi mdrun -deffnm minim1

