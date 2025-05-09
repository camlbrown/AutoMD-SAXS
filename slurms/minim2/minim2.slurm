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

cd $MINIM2_DIR && cp $MINIM1_DIR/minim1.gro . && cp $MINIM1_DIR/topol.top .

gmx_mpi grompp -f $MDP_DIR/minim2.mdp -c minim1.gro -p topol.top -o minim2.tpr -maxwarn 5 2>&1 | tee minim2_grompp.log

gmx_mpi mdrun -deffnm minim2
