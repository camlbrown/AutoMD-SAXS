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

cd $NVT_DIR && cp $MINIM2_DIR/minim2.gro . && cp $MINIM2_DIR/topol.top . && cp $MINIM2_DIR/*.itp .

gmx_mpi make_ndx -f minim2.gro -o lig.ndx << EOF
q
EOF

#output gro file with 'Other' (i.e., the ligand). 
gmx_mpi editconf -f minim2.gro -n lig.ndx -o lig.gro << EOF 
Other
EOF


gmx_mpi make_ndx -f lig.gro -o index_lig.ndx << EOF
0 & ! a H*
q
EOF

gmx_mpi genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000 << EOF
System_&_!H*
EOF

file_path='topol.top'

awk '/#include "Ligand.itp"/ {print; print ""; print "; Ligand position restraints"; print "#ifdef POSRES"; print "#include \"posre_lig.itp\""; print "#endif"; next}1' "$file_path" > temp && mv temp "$file_path"

gmx_mpi make_ndx -f minim2.gro -o index.ndx << EOF
"Protein" | "Other"
q
EOF

gmx_mpi grompp -f $MDP_DIR/lig_nvt.mdp -c minim2.gro -r minim2.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 100 2>&1 | tee nvt_grompp.log

gmx_mpi mdrun -deffnm nvt

