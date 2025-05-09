#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=02:00:00
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

cd $NPT_DIR && cp $NVT_DIR/nvt.gro . && cp $NVT_DIR/topol.top . && cp $NVT_DIR/nvt.cpt . && cp $NVT_DIR/*.itp . && cp $NVT_DIR/index.ndx . 

gmx_mpi grompp -f $MDP_DIR/lig_npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index.ndx -maxwarn 100 2>&1 | tee npt_grompp.log

gmx_mpi mdrun -deffnm npt

