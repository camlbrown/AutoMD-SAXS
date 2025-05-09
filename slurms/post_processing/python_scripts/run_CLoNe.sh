#!/bin/bash
				##################################CLONE CLUSTERING ANALYSIS###################################


#Will need to either set paths for the CLoNe scripts or copy into wk dir 
#Requires final.xtc from each repeat and final.pdb. repeat trajectories will need to be renamed to include *1 *2 *3

#make directory for combined trajectories and insert GMX cat commands here 


#cat trajectories
gmx trjcat -f *1.xtc *2.xtc *3.xtc -o combined.xtc

wait 

#align to starting structure
gmx trjconv -s $PDB2GMX_DIR/GMX.gro -f combined.xtc -o combined_aligned.xtc -fit rot+trans << EOF
3
1
EOF

wait

#Run CLoNe clustering
wait
mkdir pdc1 && cd pdc1
python ../run_structural.py -traj ../combined_aligned.xtc -topo ../GMX.gro -pdc 1 -at_sel "name CA" -pca 2 << EOF 
3
EOF

wait

cd .. && mkdir pdc3 && cd pdc3 
python ../run_structural.py -traj ../combined_aligned.xtc -topo ../GMX.gro -pdc 3 -at_sel "name CA" -pca 2 << EOF
3
EOF

wait

cd .. && mkdir pdc5 && cd pdc5
python ../run_structural.py -traj ../combined_aligned.xtc -topo ../GMX.gro -pdc 5 -at_sel "name CA" -pca 2 << EOF
3
EOF

wait

cd .. && mkdir pdc7 && cd pdc7 
python ../run_structural.py -traj ../combined_aligned.xtc -topo ../GMX.gro -pdc 7 -at_sel "name CA" -pca 2 << EOF
3
EOF


