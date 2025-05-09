#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=2
#SBATCH --time=20:00:00
#SBATCH --partition=batch

# Check if SIMULATION_DIR is set
if [ -z "$SIMULATION_DIR" ]; then
    echo "ERROR: SIMULATION_DIR is not recognised by post processing script. link between scripts needs fixing."
    exit 1
else
    echo "SIMULATION_DIR is set to: $SIMULATION_DIR" 
fi

config_file="$SIMULATION_DIR/configurations.txt"

if [ ! -f "$config_file" ]; then
    echo "ERROR: configurations.txt file not found in $SIMULATION_DIR."
    exit 1
fi

source "$config_file"

export GMXLIB="$BASE_DIR/ff_files"


##############################
#REPEAT1 TRAJECTORY PROCESSING
##############################

load_gmx

cd $PROCESSED_R1 && cp $REPEAT_DIR1/*.gro . && cp $REPEAT_DIR1/*.top . && cp $REPEAT_DIR1/*.tpr . 

gmx_mpi editconf -f *.gro -o model.pdb -c


gmx_mpi make_ndx -f model.pdb -o nosolv.ndx << EOF 
"Protein" | "Other"
q
EOF

#remove solvent from gro, tpr and xtc 
gmx_mpi editconf -f model.pdb -n *ndx -o nosolv.gro << EOF
Protein_Other
EOF

gmx_mpi convert-tpr -s *.tpr -n *.ndx -o nosolv.tpr << EOF
Protein_Other
EOF

gmx_mpi trjconv -f ../*.xtc -s nosolv.tpr -n *ndx -o nosolv.xtc <<EOF
Protein_Other
EOF

#nojump
gmx_mpi trjconv -f nosolv.xtc -s nosolv.tpr  -pbc nojump -o nojump.xtc << EOF
0
EOF

gmx_mpi trjconv -f nojump.xtc -s nosolv.gro -o final.xtc -fit rot+trans << EOF
0
0
EOF

mkdir final_pdb && cd final_pdb && cp $PDB2GMX_DIR/Complex.pdb .
mv Complex.pdb final.pdb 

cp final.pdb ../

cd $R1_EXTRACT_FRAMES && cp ../final_pdb/final.pdb . && cp ../final.xtc . && cp ../nosolv.tpr .

wait 

#extract frames
gmx_mpi trjconv -f final.xtc -s nosolv.tpr -o structure_.pdb -sep -skip 2 << EOF
1
EOF

##########################
#NON-SAXS MD ANALYSIS REP1
##########################

cd $PROCESSED_R1

# C-alpha RMSD
gmx_mpi rms -s nosolv.tpr -f final.xtc -o rmsd_backbone.xvg -tu ns << EOF
3
3
EOF

# Gyration
gmx_mpi gyrate -s nosolv.tpr -f final.xtc -o gyrate.xvg << EOF
System
EOF

# SASA
gmx_mpi sasa -s nosolv.tpr -f final.xtc -o sasa_total.xvg -surface Protein -output Protein

# System hbonds 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_system.xvg << EOF
System
System
EOF

# Total energy
gmx_mpi energy -f ../md.edr -o energy.xvg << EOF
Total-Energy

EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/standard_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir standard_analysis
mv *svg *pdf standard_analysis
mv xvg_files standard_analysis

load_gmx

################
#Ligand analysis
################

gmx_mpi make_ndx -f final.pdb -o ligand.ndx << EOF
q
EOF

#rmsd ligand 
gmx_mpi rms -s nosolv.tpr -f final.xtc -o rmsd_ligand.xvg -n ligand.ndx << EOF
Protein 
Other
EOF

#rmsf ligand 
gmx_mpi rmsf -s nosolv.tpr -f final.xtc -o rmsf_ligand.xvg -n ligand.ndx << EOF
Other
EOF

#hbonds between prot and lig 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_lig.xvg << EOF
Protein
Other
EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/prot-lig_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir prot-lig_analysis
mv *svg *pdf prot-lig_analysis
mv xvg_files prot-lig_analysis

load_gmx

##############################
#REPEAT2 TRAJECTORY PROCESSING
##############################

#create nosolvent index file
cd $PROCESSED_R2 && cp $REPEAT_DIR2/*.gro . && cp $REPEAT_DIR2/*.top . && cp $REPEAT_DIR2/*.tpr . 

gmx_mpi editconf -f *.gro -o model.pdb -c

gmx_mpi make_ndx -f model.pdb -o nosolv.ndx << EOF 
"Protein" | "Other"
q
EOF

#remove solvent from gro, tpr and xtc 
gmx_mpi editconf -f model.pdb -n *ndx -o nosolv.gro << EOF
Protein_Other
EOF

gmx_mpi convert-tpr -s *.tpr -n *.ndx -o nosolv.tpr << EOF
Protein_Other
EOF

gmx_mpi trjconv -f ../*.xtc -s nosolv.tpr -n *ndx -o nosolv.xtc <<EOF
Protein_Other
EOF

#nojump
gmx_mpi trjconv -f nosolv.xtc -s nosolv.tpr  -pbc nojump -o nojump.xtc << EOF
0
EOF

gmx_mpi trjconv -f nojump.xtc -s nosolv.gro -o final.xtc -fit rot+trans << EOF
0
0
EOF

mkdir final_pdb && cd final_pdb && cp $PDB2GMX_DIR/Complex.pdb .
mv Complex.pdb final.pdb 

cp final.pdb ../

cd $R2_EXTRACT_FRAMES && cp ../final_pdb/final.pdb . && cp ../final.xtc . && cp ../nosolv.tpr .

#extract frames
gmx_mpi trjconv -f final.xtc -s nosolv.tpr -o structure_.pdb -sep -skip 2 << EOF
1
EOF

##########################
#NON-SAXS MD ANALYSIS REP2
##########################

cd $PROCESSED_R2

# C-alpha RMSD
gmx_mpi rms -s nosolv.tpr -f final.xtc -o rmsd_backbone.xvg -tu ns << EOF
3
3
EOF

# Gyration
gmx_mpi gyrate -s nosolv.tpr -f final.xtc -o gyrate.xvg << EOF
System
EOF

# SASA
gmx_mpi sasa -s nosolv.tpr -f final.xtc -o sasa_total.xvg -surface Protein -output Protein

# System hbonds 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_system.xvg << EOF
System
System
EOF

# Total energy
gmx_mpi energy -f ../md.edr -o energy.xvg << EOF
Total-Energy

EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/standard_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir standard_analysis
mv *svg *pdf standard_analysis
mv xvg_files standard_analysis

load_gmx

################
#Ligand analysis
################

gmx_mpi make_ndx -f final.pdb -o ligand.ndx << EOF
q
EOF

#rmsd ligand 
gmx_mpi rms -s nosolv.tpr -f final.xtc -o rmsd_ligand.xvg -n ligand.ndx << EOF
Protein 
Other
EOF

#rmsf ligand 
gmx_mpi rmsf -s nosolv.tpr -f final.xtc -o rmsf_ligand.xvg -n ligand.ndx << EOF
Other
EOF

#hbonds between prot and lig 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_lig.xvg << EOF
Protein
Other
EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/prot-lig_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir prot-lig_analysis
mv *svg *pdf prot-lig_analysis
mv xvg_files prot-lig_analysis

load_gmx

##############################
#REPEAT3 TRAJECTORY PROCESSING
##############################

cd $PROCESSED_R3 && cp $REPEAT_DIR3/*.gro . && cp $REPEAT_DIR3/*.top . && cp $REPEAT_DIR3/*.tpr . 

gmx_mpi editconf -f *.gro -o model.pdb -c

gmx_mpi make_ndx -f model.pdb -o nosolv.ndx << EOF 
"Protein" | "Other"
q
EOF

#remove solvent from gro, tpr and xtc 
gmx_mpi editconf -f model.pdb -n *ndx -o nosolv.gro << EOF
Protein_Other
EOF

gmx_mpi convert-tpr -s *.tpr -n *.ndx -o nosolv.tpr << EOF
Protein_Other
EOF

#mv $REPEAT_DIR1/*xtc .
gmx_mpi trjconv -f ../*.xtc -s nosolv.tpr -n *ndx -o nosolv.xtc <<EOF
Protein_Other
EOF

#nojump
gmx_mpi trjconv -f nosolv.xtc -s nosolv.tpr  -pbc nojump -o nojump.xtc << EOF
0
EOF

gmx_mpi trjconv -f nojump.xtc -s nosolv.gro -o final.xtc -fit rot+trans << EOF
0
0
EOF

mkdir final_pdb && cd final_pdb && cp $PDB2GMX_DIR/Complex.pdb .
mv Complex.pdb final.pdb 

cp final.pdb ../

cd $R3_EXTRACT_FRAMES && cp ../final_pdb/final.pdb . && cp ../final.xtc . && cp ../nosolv.tpr .

#extract frames
gmx_mpi trjconv -f final.xtc -s nosolv.tpr -o structure_.pdb -sep -skip 2 << EOF
1
EOF

##########################
#NON-SAXS MD ANALYSIS REP3
##########################

cd $PROCESSED_R3

# C-alpha RMSD
gmx_mpi rms -s nosolv.tpr -f final.xtc  -o rmsd_backbone.xvg -tu ns << EOF
3
3
EOF

# Gyration
gmx_mpi gyrate -s nosolv.tpr -f final.xtc -o gyrate.xvg << EOF
System
EOF

# SASA
gmx_mpi sasa -s nosolv.tpr -f final.xtc -o sasa_total.xvg -surface Protein -output Protein

# System hbonds 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_system.xvg << EOF
System
System
EOF

# Total energy
gmx_mpi energy -f ../md.edr -o energy.xvg << EOF
Total-Energy

EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/standard_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir standard_analysis
mv *svg *pdf standard_analysis
mv xvg_files standard_analysis

################
#Ligand analysis
################

gmx_mpi make_ndx -f final.pdb -o ligand.ndx << EOF
q
EOF

#rmsd ligand 
gmx_mpi rms -s nosolv.tpr -f final.xtc -o rmsd_ligand.xvg -n ligand.ndx << EOF
Protein 
Other
EOF

#rmsf ligand 
gmx_mpi rmsf -s nosolv.tpr -f final.xtc -o rmsf_ligand.xvg -n ligand.ndx << EOF
Other
EOF

#hbonds between prot and lig 
gmx_mpi hbond -s nosolv.tpr -f final.xtc -num hb_lig.xvg << EOF
Protein
Other
EOF

unload_gmx

$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/prot-lig_md_analysis.py

mkdir xvg_files
mv *xvg xvg_files
mkdir prot-lig_analysis
mv *svg *pdf prot-lig_analysis
mv xvg_files prot-lig_analysis

load_gmx

#################
#CLoNe clustering 
#################

cd $PRODUCTION_DIR 
mkdir combined_traj 
cd combined_traj

cp $SLURM_DIR/post_processing/python_scripts/* .

load_gmx

#cat trajectories
echo 'Combining trajectory repeats for clustering analysis'
gmx_mpi trjcat -f $PROCESSED_R1/final.xtc $PROCESSED_R2/final.xtc $PROCESSED_R3/final.xtc -o combined.xtc -cat

wait 

#align to starting structure
gmx_mpi trjconv -s $PROCESSED_R1/final_pdb/final.pdb -f combined.xtc -o combined_aligned.xtc -skip 2 -fit rot+trans << EOF
3
0
EOF

wait

#Run CLoNe clustering
wait
echo 'Running CLoNE clustering'

mkdir clustering && cd clustering
unload_gmx

mkdir pdc1 && cd pdc1
$PYTHON_CMD $PRODUCTION_DIR/combined_traj/run_structural.py -traj $PRODUCTION_DIR/combined_traj/combined_aligned.xtc -topo $PROCESSED_R1/final_pdb/final.pdb -pdc 1 -at_sel "name CA" -pca 2 << EOF
3
EOF

wait

cd .. && mkdir pdc3 && cd pdc3 
$PYTHON_CMD $PRODUCTION_DIR/combined_traj/run_structural.py -traj $PRODUCTION_DIR/combined_traj/combined_aligned.xtc -topo $PROCESSED_R1/final_pdb/final.pdb -pdc 3 -at_sel "name CA" -pca 2 << EOF
3
EOF

wait

cd .. && mkdir pdc5 && cd pdc5
$PYTHON_CMD $PRODUCTION_DIR/combined_traj/run_structural.py -traj $PRODUCTION_DIR/combined_traj/combined_aligned.xtc -topo $PROCESSED_R1/final_pdb/final.pdb -pdc 5 -at_sel "name CA" -pca 2 << EOF
3
EOF

wait

cd .. && mkdir pdc7 && cd pdc7 
$PYTHON_CMD $PRODUCTION_DIR/combined_traj/run_structural.py -traj $PRODUCTION_DIR/combined_traj/combined_aligned.xtc -topo $PROCESSED_R1/final_pdb/final.pdb -pdc 7 -at_sel "name CA" -pca 2 << EOF
3
EOF

#Invoke SAXS analysis if SAXS file present 

if [[ "$SAXS_FILE" != "None" ]]; then
  echo "Running SAXS analysis..."
  sh $SLURM_DIR/post_processing/saxs_analysis.sh "$SIMULATION_DIR"
else
  echo "No SAXS data provided; skipping SAXS analysis."
fi
