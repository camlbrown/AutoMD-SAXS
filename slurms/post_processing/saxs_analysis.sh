#########SAXS ANALYSIS SCRIPT##########

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


# RENAME RESIDUES FOR CRYSOL

cd $R1_EXTRACT_FRAMES
#process extracted structures for SAXS analysis (CRYSOL naming conventions)
# Remove caps
echo 'Removing caps'
sed -i '/ACE/d' structure_*
sed -i '/NME/d' structure_*

# Check the FORCE_FIELD variable and run the corresponding block
if [[ "$FORCE_FIELD" == "amber14sb" ]]; then
    # amber block
    echo 'Renaming non-standard residues for SAXS analysis'
    echo 'Renaming HIE/HID/HIP'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*

elif [[ "$FORCE_FIELD" == "charmm36m" ]]; then
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
else
    echo "Unknown FORCE_FIELD: $FORCE_FIELD. Running both amber and charmm blocks."
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis (charmm block)'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
fi
wait

cd $R2_EXTRACT_FRAMES

#process extracted structures for SAXS analysis (CRYSOL naming conventions)
# Remove caps
echo 'Removing caps'
sed -i '/ACE/d' structure_*
sed -i '/NME/d' structure_*

# Check the FORCE_FIELD variable and run the corresponding block
if [[ "$FORCE_FIELD" == "amber14sb" ]]; then
    # amber block
    echo 'Renaming non-standard residues for SAXS analysis'
    echo 'Renaming HIE/HID/HIP'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*

elif [[ "$FORCE_FIELD" == "charmm36m" ]]; then
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
else
    echo "Unknown FORCE_FIELD: $FORCE_FIELD. Running both amber and charmm blocks."
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis (charmm block)'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
fi
wait


cd $R3_EXTRACT_FRAMES

#process extracted structures for SAXS analysis (CRYSOL naming conventions)
# Remove caps
echo 'Removing caps'
sed -i '/ACE/d' structure_*
sed -i '/NME/d' structure_*

# Check the FORCE_FIELD variable and run the corresponding block
if [[ "$FORCE_FIELD" == "amber14sb" ]]; then
    # amber block
    echo 'Renaming non-standard residues for SAXS analysis'
    echo 'Renaming HIE/HID/HIP'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*

elif [[ "$FORCE_FIELD" == "charmm36m" ]]; then
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
else
    echo "Unknown FORCE_FIELD: $FORCE_FIELD. Running both amber and charmm blocks."
    echo 'Renaming non-standard residues for SAXS analysis'
    sed -i 's/HIE/HIS/g' structure_*
    sed -i 's/HID/HIS/g' structure_*
    sed -i 's/HIP/HIS/g' structure_*
    sed -i 's/ASH/ASP/g' structure_*
    sed -i 's/GLH/GLU/g' structure_*
    sed -i 's/CYX/CYS/g' structure_*
    # charmm block
    echo 'Renaming non-standard residues for SAXS analysis (charmm block)'
    sed -i 's/HSE/HIS/g' structure_*
    sed -i 's/HSD/HIS/g' structure_*
    sed -i 's/HSP/HIS/g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/ASPP/ASP /g' structure_*
    sed -i 's/GLUP/GLU /g' structure_*
    sed -i 's/CD  ILE/CD1 ILE/g' structure_*
    sed -i '/OT1\|OT2/d' structure_*
fi
wait

# SAXS CALCULATIONS 

# shanum
shanum_output="$(shanum "$SAXS_DIR"/*.dat 2>&1)"

# Extract the first line containing "Smax="
SMAX=$(printf '%s\n' "$shanum_output" \
    | grep 'Smax=' \
    | head -1 \
    | sed -E 's/.*Smax=[[:space:]]*//')

SMAX=$(printf "%.2f" "$SMAX")
export SMAX
echo "Detected Smax = $SMAX"


cd $R1_EXTRACT_FRAMES 
crysol $SAXS_DIR/*dat structure* -lm 30 -cst -sm $SMAX 
wait
mkdir crysol_r1
mkdir fit_log 
mv *fit *log fit_log
mv fit_log crysol_r1 
awk '/Model:/ { found=1 } found { print }' crysol_summary.txt > temp && mv temp crysol_summary.txt
mv crysol_summary.txt crysol_r1
cd crysol_r1
$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/chi2vsRg.py 

cd $R2_EXTRACT_FRAMES 
crysol $SAXS_DIR/*dat structure* -lm 30 -cst -sm $SMAX 
wait
mkdir crysol_r2
mkdir fit_log 
mv *fit *log fit_log
mv fit_log crysol_r2
awk '/Model:/ { found=1 } found { print }' crysol_summary.txt > temp && mv temp crysol_summary.txt
mv crysol_summary.txt crysol_r2
cd crysol_r2
$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/chi2vsRg.py 

cd $R3_EXTRACT_FRAMES 
crysol $SAXS_DIR/*dat structure* -lm 30 -cst -sm $SMAX   
wait
mkdir crysol_r3
mkdir fit_log 
mv *fit *log fit_log
mv fit_log crysol_r3
awk '/Model:/ { found=1 } found { print }' crysol_summary.txt > temp && mv temp crysol_summary.txt
mv crysol_summary.txt crysol_r3
cd crysol_r3
$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/chi2vsRg.py

# SAXS-SCORED PCA

cd $PRODUCTION_DIR/combined_traj  

mkdir extract_frames && cd extract_frames 

load_gmx

gmx_mpi trjconv -f $PRODUCTION_DIR/combined_traj/combined_aligned.xtc -s $PROCESSED_R1/final_pdb/final.pdb -o structure_.pdb -sep << EOF
0
EOF

wait

cp $SAXS_DIR/*dat . && cp $PROCESSED_R1/final_pdb/final.pdb .

unload_gmx

crysol $SAXS_DIR/*dat structure* -lm 30 -cst -sm $SMAX 

wait

awk '/Model:/ { found=1 } found { print }' crysol_summary.txt > temp && mv temp crysol_summary.txt

mkdir fit_log
mv *fit *log fit_log
mkdir crysol
mv fit_log crysol
cp crysol_summary.txt crysol
cd crysol
$PYTHON_CMD $SLURM_DIR/post_processing/python_scripts/PCA_viridis.py 

cd ../

cp $SAXS_DIR/*dat .

# GAJOE 
# using smaxs from shanum, 30 spherical harms, else default
gajoe -p << EOF
1
$SAXS_FILE
.
30
$SMAX 
51
1000
50
no
20
5
yes
yes
100
EOF

wait

mv junXXX.eom Size_listXXX.txt RanchXXX.log GA001

cd GA001/curve_1

~/.conda/envs/md/bin/python $SLURM_DIR/post_processing/python_scripts/plot_rg_ensemble.py

~/.conda/envs/md/bin/python $SLURM_DIR/post_processing/python_scripts/ensemble_scattering.py
