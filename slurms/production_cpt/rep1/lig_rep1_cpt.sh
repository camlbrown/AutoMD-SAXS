#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=01:10:00
#SBATCH --partition=batch

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

cd "$REPEAT_DIR1"

mpirun gmx_mpi mdrun -s md.tpr -cpi md.cpt -deffnm md -maxh 1 
wait

sleep 100

if [[ -n "$EMAIL_ADDR" ]]; then
  MAIL_FLAGS="--mail-type=ALL --mail-user=$EMAIL_ADDR"
else
  MAIL_FLAGS=""
fi


# Check if .gro file exists
if [ -e *.gro ]; then
    if grep -q "post_processing initiated" $SIMULATION_DIR/monitor.txt; then
        echo "A post-processing job has already been initiated."
    else
        echo "repeat1 completed" >> $SIMULATION_DIR/monitor.txt
        if grep -q "repeat1 completed" $SIMULATION_DIR/monitor.txt && \
           grep -q "repeat2 completed" $SIMULATION_DIR/monitor.txt && \
           grep -q "repeat3 completed" $SIMULATION_DIR/monitor.txt; then
            echo "All repeats have finished. Initiating post-processing..."
            JOBID_POST=$(sbatch --parsable $MAIL_FLAGS -J post --export=ALL,SIMULATION_DIR="$SIMULATION_DIR" "$SLURM_DIR/post_processing/prot-lig_post_processing.sh" "$SIMULATION_DIR" 2>&1)
            echo "JOBID_POST=${JOBID_POST}" >> $SIMULATION_DIR/post_processing.log
            echo "post_processing initiated" >> $SIMULATION_DIR/monitor.txt
        fi
    fi
else
    # Resubmit script if .gro file is not found
    JOBID_MDR1_cpt=$(sbatch --parsable -J lr1_cpt --export=ALL "$SLURM_DIR/production_cpt/rep1/lig_rep1_cpt.sh" "$SIMULATION_DIR")
    echo "JOBID_MDR1_cpt=${JOBID_MDR1_cpt}"
fi
