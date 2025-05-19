#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=batch

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

cd $NVT_DIR && cp $MINIM2_DIR/minim2.gro . && cp $MINIM2_DIR/topol.top . 

cp $PDB2GMX_dir/posre.itp .

gmx_mpi grompp -f $MDP_DIR/nvt.mdp -c minim2.gro -r minim2.gro -p topol.top -o nvt.tpr 2>&1 | tee nvt_grompp.log

gmx_mpi mdrun -deffnm nvt

