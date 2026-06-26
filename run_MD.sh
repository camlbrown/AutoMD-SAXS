#!/bin/bash

#-----------------MD RUN SCRIPT-----------------

# General housekeeping
set -euo pipefail


# Parse arguments. A new optional --dry-run flag delegates validation and
# command planning to the Python package (automd_saxs) without running GROMACS,
# ATSAS, or Slurm. Everything else is unchanged.
DRY_RUN=0
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    -h|-help|--help)
      cat <<EOM
Automated MD with GROMACS

Usage: bash run_MD.sh <simulation_directory> [--dry-run]

Run only after setting up your simulation with simulation_setup.sh.
  --dry-run   Validate the configuration and print the planned commands
              (GROMACS setup, Slurm DAG, directories) without executing
              anything. Requires the automd_saxs Python package on PYTHONPATH.
EOM
      exit 0
      ;;
    --dry-run) DRY_RUN=1 ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done

if [ "${#POSITIONAL[@]}" -ne 1 ]; then
    echo "Usage: $0 <simulation_directory> [--dry-run]"
    exit 1
fi
set -- "${POSITIONAL[@]}"

SIMULATION_DIR="$1"

# Source configurations.txt in simulation directory
config_file="$SIMULATION_DIR/configurations.txt"
if [ -f "$config_file" ]; then
    source "$config_file"
else
    echo "Error: configurations.txt file not found in $SIMULATION_DIR."
    exit 1
fi

GMXLIB="$BASE_DIR/ff_files"
export GMXLIB
FF_NAME=${FORCE_FIELD##*/}

# Dry-run: delegate validation and command planning to the Python package and
# exit before touching GROMACS/Slurm. Lets the workflow be inspected on machines
# without the scientific binaries installed.
if [ "$DRY_RUN" -eq 1 ]; then
  DRY_PDB="$SIMULATION_DIR/pdb2gmx/GMX.pdb"
  PDB_ARG=()
  [ -f "$DRY_PDB" ] && PDB_ARG=(--pdb "$DRY_PDB")
  PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
    -m automd_saxs plan \
    --config "$config_file" \
    --work-dir "$BASE_DIR" \
    --slurm-dir "$SLURM_DIR" \
    --mdp-dir "$MDP_DIR" \
    "${PDB_ARG[@]}"
  exit $?
fi

# Conditional setup

# Check if DISULFIDE is set to 'y' or 'yes'
if [[ "$DISULFIDE" == "y" || "$DISULFIDE" == "yes" ]]; then
    ss_flag="-ss"
else
    ss_flag=""
fi

# Choose the production .mdp file based on SYSTEM
if [ "$SYSTEM" = "Protein-ligand" ]; then
    TIMESTEP="$MDP_DIR/lig_md.mdp"
elif [ "$SYSTEM" = "Protein" ]; then
    TIMESTEP="$MDP_DIR/md.mdp"
else
    echo "Error: Unknown SYSTEM type: $SYSTEM" >&2
    exit 1
fi

# Set the production step count (and dt) from the configuration. Rendered with
# the Python package instead of an in-place `sed -i`, which previously mutated
# the git-tracked template on every run. $MDP_DIR is now a per-job copy (created
# by simulation_setup.sh), so rendering in place here is safe and reproducible.
PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
  -m automd_saxs render-mdp \
  --template "$TIMESTEP" \
  --out "$TIMESTEP" \
  --config "$config_file"
echo "Rendered production parameters into $(basename "$TIMESTEP") (nsteps from SIMULATION_TIME=$SIMULATION_TIME ns)"


if [[ "$SYSTEM" == "Protein-ligand" ]]; then

    # copy and split out ligand/protein
    cd "$LIGAND_SETUP"
    mv GMX.pdb complex.pdb
    $PYTHON_CMD \
      "$SLURM_DIR/ligand_setup/separate_multi_lig.py"

    # prepare protein
    cp protein.pdb ProteinAmber.pdb
    load_gmx
    gmx_mpi pdb2gmx -ff amber14sb -f ProteinAmber.pdb -o Protein_pdb2gmx.pdb -p Protein.top -ter -water spce -ignh $ss_flag
    unload_gmx

    cp Protein_pdb2gmx.pdb Complex.pdb
    sed -i '/ENDMDL/d' Complex.pdb

    get_chain_identifier() {
        ligand_filename=$1
        chain_identifier=$(echo "$ligand_filename" | sed -n 's/.*_\([A-Za-z]\).pdb/\1/p')
        echo "$chain_identifier"
    }

# Iterate over ligand files in numerical order
    for ligand_file in $(ls -v $LIGAND_SETUP/ligands/ligand_*.pdb); do
        ligand_name=$(basename "$ligand_file")
        ligand_name="${ligand_name%.*}" 
        cp "$ligand_file" "${ligand_name}_H.pdb"
        reduce "$ligand_file" > "${ligand_name}_H.pdb" # Add hydrogens to ligand
        wait  
        pdb4amber -i "${ligand_name}_H.pdb" -o "${ligand_name}.pdb" # clean 
        wait
        acpype -a gaff2 -i "${ligand_name}.pdb" -b "$ligand_name" # parameterise 
        wait
 
        ligand_acpype_dir="${ligand_name}.acpype"
        python $SLURM_DIR/ligand_setup/replace_chain_identifier.py "${ligand_file}" "$LIGAND_SETUP/${ligand_acpype_dir}/${ligand_name}_NEW.pdb"
        grep -h ATOM "$LIGAND_SETUP/${ligand_acpype_dir}/${ligand_name}_NEW.pdb" >> Complex.pdb
        cp "${ligand_acpype_dir}/${ligand_name}_GMX.itp" "${ligand_name}.itp" 
        echo "copy protein.top to complex_ligname.top"
        cp Protein.top Complex_${ligand_name}.top
        insert_text="#include \"${ligand_name}.itp\""
        input_file="Protein.top"
        sed -i '/#include "amber.*\/forcefield.itp"/ {
            a\
            \
'"$insert_text"'
        }' "$input_file"

        cat "$input_file" | sed '/forcefield.itp\"/a\' >| "${input_file}_2.top"
        echo "$ligand_name         1" >> "${input_file}_2.top"
        mv "${input_file}_2.top" "$input_file"
    done
    echo "moving protein.top to complex.top"

    mv Protein.top topol.top
    cp topol.top "$PDB2GMX_DIR"
    cp Complex.pdb "$PDB2GMX_DIR"

    shopt -s nullglob

    itps=( "$LIGAND_SETUP"/*.itp )
    if (( ${#itps[@]} )); then
        cp "$LIGAND_SETUP"/*.itp "$MINIM1_DIR"
        cp "$LIGAND_SETUP"/*.itp "$MINIM2_DIR"
        cp "$LIGAND_SETUP"/*.itp "$NVT_DIR"
        cp "$LIGAND_SETUP"/*.itp "$NPT_DIR"
        cp "$LIGAND_SETUP"/*.itp "$PRODUCTION_DIR"
        cp "$LIGAND_SETUP"/*.itp "$GENION_DIR"
        cp "$LIGAND_SETUP"/*.itp "$SOLVATE_DIR"
    fi

    cd "$PDB2GMX_DIR"

    pdb_tidy Complex.pdb \
      >> Complex_tidy.pdb
    $PYTHON_CMD \
      "$SLURM_DIR/ligand_setup/ATOM_to_HETATM.py"
    load_gmx
    gmx_mpi editconf -f GMX.pdb -o GMX.gro 
    unload_gmx



elif [[ "$SYSTEM" == "Protein" ]]; then

    load_gmx
    cd $PDB2GMX_DIR
    gmx_mpi pdb2gmx -f GMX.pdb -o GMX.gro -p topol.top -chainsep ter -ff $FF_NAME -water tip3p -ter -merge all $ss_flag

    shopt -s nullglob

    itps=( "$PDB2GMX_DIR"/*.itp )
    if (( ${#itps[@]} )); then
        cp "$PDB2GMX_DIR"/*.itp "$MINIM1_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$MINIM2_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$NVT_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$NPT_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$PRODUCTION_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$GENION_DIR"
        cp "$PDB2GMX_DIR"/*.itp "$SOLVATE_DIR"
    fi

fi

#BOX SETUP

cd $SOLVATE_DIR && cp $PDB2GMX_DIR/GMX.gro . && cp $PDB2GMX_DIR/*top . 

if [[ "$SAXS_FILE" != "None" ]]; then
  cp $SAXS_DIR/*dat .
else
  echo "No SAXS file, using model Dmax as buffer"
fi

DMAX_PDB_FILE="$SIMULATION_DIR/pdb2gmx/GMX.pdb"

# Compute model Dmax and solvation box padding via the Python package. This
# replaces an inline numpy/scipy heredoc plus `bc` arithmetic (both untested and
# duplicated). The subcommand prints MODEL_DMAX / EFFECTIVE_DMAX / BOX_PADDING as
# numeric shell assignments; "$DMAX" may be a number or the literal "Model".
dmax_out=$(PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON_CMD:-python3}" \
  -m automd_saxs dmax --pdb "$DMAX_PDB_FILE" --experimental "$DMAX") || {
    echo "Error: Could not determine box padding from $DMAX_PDB_FILE" >&2
    exit 1
}
eval "$dmax_out"
DMAX="$EFFECTIVE_DMAX"

echo "Model Dmax:        $MODEL_DMAX nm"
echo "Experimental Dmax: $DMAX nm"
echo "Using box padding: $BOX_PADDING nm"


load_gmx

# Box selection based on protein shape
if [ "$BOX_SHAPE" == "auto" ]; then
    gmx_mpi gyrate -f *.gro -s *.tpr -o gyrate.xvg
    RGX=$(awk '$1 ~ /^[0-9]/ {print $2}' gyrate.xvg | sort -nr | head -1)  # max rg 

    if (( $(echo "$RGX < 2.0" | bc -l) )); then
        BOX_SHAPE="dodecahedron"  # for compact proteins
    else
        BOX_SHAPE="triclinic"  # for elongated shapes
    fi
fi
echo "$BOX_SHAPE"

case $BOX_SHAPE in
    dodecahedron)
        gmx_mpi editconf -f GMX.gro -o 1.gro -bt dodecahedron -d $BOX_PADDING 2>&1 | tee editconf.log
        ;;
    octahedron)
        gmx_mpi editconf -f GMX.gro -o 1.gro -bt octahedron -d $BOX_PADDING 2>&1 | tee editconf.log
        ;;
    triclinic|rectangular)
        gmx_mpi editconf -f GMX.gro -o 1.gro -bt triclinic -d $BOX_PADDING 2>&1 | tee editconf.log
        ;;
    cubic)
        gmx_mpi editconf -f GMX.gro -o 1.gro -bt cubic -d $BOX_PADDING 2>&1 | tee editconf.log
        ;;
    *)
        echo "Invalid BOX_SHAPE setting."
        exit 1
        ;;
esac

read x_box y_box z_box <<< $(awk 'END{print $1, $2, $3}' 1.gro)

# Centre protein
gmx_mpi editconf -f 1.gro -o centered.gro -c 2>&1 | tee editconf.log # centering with -c will become apparant after trjconv with tpr input 

# Solvate
gmx_mpi solvate -cp centered.gro -cs spc216.gro -o solvate.gro -p *.top 2>&1 | tee solvate.log

cd $GENION_DIR && cp $SOLVATE_DIR/solvate.gro . && cp $SOLVATE_DIR/*.top .

# conditional prot-lig itp handling
if [ "$SYSTEM" = "Protein-ligand" ]; then

    # 1) Iterate over n ligand itp files
    for ligand_itp_file in ligand_*.itp; do
        # Extract the ligand name from the itp file
        ligand_name=$(basename "$ligand_itp_file" .itp)

        # Find the line number where force field inclusion starts in top
        force_field_line=$(grep -n '#include ".*forcefield.itp"' topol.top | cut -d: -f1)

        # Find the line numbers for the atomtypes section in the ligand itp
        start_line=$(grep -n '\[ atomtypes \]' "$ligand_itp_file" | cut -d: -f1)
        end_line=$(grep -n '\[ moleculetype \]' "$ligand_itp_file" | cut -d: -f1)
        end_line=$((end_line - 1))  # exclude the [ moleculetype ] line

        # Extract and remove the atomtypes section
        sed -n "${start_line},${end_line}p" "$ligand_itp_file" > "${ligand_name}_atomtypes.txt"
        sed -i "${start_line},${end_line}d" "$ligand_itp_file"

        # Append the extracted atomtypes into top
        sed -i "${force_field_line}a\\
#include \"${ligand_name}_atomtypes.txt\" " topol.top
    done
    echo "itp stuff done, now copy"
    # 2) Copy all .txt and .itp files into each run directory
    for dir in "$MINIM1_DIR" "$MINIM2_DIR" "$NVT_DIR" "$NPT_DIR" "$PRODUCTION_DIR"; do
        cp *.txt "$dir"
        cp *.itp "$dir"
    done

fi

# Neutralisation

gmx_mpi grompp -f $MDP_DIR/ions.mdp -c *.gro -p *.top -o ions.tpr -maxwarn 100 

# Use the user-selected ionic concentration. The legacy script hardcoded
# -conc 0.15, silently ignoring IONIC_CONCENTRATION from configurations.txt.
echo "SOL" | gmx_mpi genion -s ions.tpr -o *.gro -p *.top -pname Na -nname Cl -neutral -conc "${IONIC_CONCENTRATION:-0.15}" 2>&1 | tee genion_grompp.log

echo "==> Submitting jobs"
# Energy minimisation
# Using steepest descent
JOBID_MINIM1=$(sbatch --parsable -J min1 $PARTITION --export=ALL "$SLURM_DIR/minim1/minim1.slurm" "$SIMULATION_DIR")
echo "JOBID_MINIM1=${JOBID_MINIM1}"

# Using conjugate gradient
JOBID_MINIM2=$(sbatch --parsable -J min2 $PARTITION --export=ALL --dependency=afterok:${JOBID_MINIM1} "$SLURM_DIR/minim2/minim2.slurm" "$SIMULATION_DIR")
echo "JOBID_MINIM2=${JOBID_MINIM2}"


# Equilibration and production

if [ "$SYSTEM" = "Protein-ligand" ]; then

  JOBID_NVT=$(sbatch --parsable \
    -J nvt \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_MINIM2} \
    "$SLURM_DIR/nvt/lig_nvt.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_NVT=${JOBID_NVT}"

  JOBID_NPT=$(sbatch --parsable \
    -J npt \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_NVT} \
    "$SLURM_DIR/npt/lig_npt.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_NPT=${JOBID_NPT}"

  JOBID_SETUP=$(sbatch --parsable \
    -J setup \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_NPT} \
    "$SLURM_DIR/production/md_setup/lig_md_setup.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_SETUP=${JOBID_SETUP}"

  # three replicates
  JOBID_MDR1=$(sbatch --parsable \
    -J rep1 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep1/lig_rep1.slurm" \
    "$SIMULATION_DIR")
  JOBID_MDR2=$(sbatch --parsable \
    -J rep2 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep2/lig_rep2.slurm" \
    "$SIMULATION_DIR")
  JOBID_MDR3=$(sbatch --parsable \
    -J rep3 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep3/lig_rep3.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_MDR1=${JOBID_MDR1}"
  echo "JOBID_MDR2=${JOBID_MDR2}"
  echo "JOBID_MDR3=${JOBID_MDR3}"

elif [ "$SYSTEM" = "Protein" ]; then

  JOBID_NVT=$(sbatch --parsable \
    -J nvt \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_MINIM2} \
    "$SLURM_DIR/nvt/nvt.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_NVT=${JOBID_NVT}"

  JOBID_NPT=$(sbatch --parsable \
    -J npt \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_NVT} \
    "$SLURM_DIR/npt/npt.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_NPT=${JOBID_NPT}"

  JOBID_SETUP=$(sbatch --parsable \
    -J setup \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_NPT} \
    "$SLURM_DIR/production/md_setup/md_setup.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_SETUP=${JOBID_SETUP}"

  JOBID_MDR1=$(sbatch --parsable \
    -J rep1 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep1/rep1.slurm" \
    "$SIMULATION_DIR")
  JOBID_MDR2=$(sbatch --parsable \
    -J rep2 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep2/rep2.slurm" \
    "$SIMULATION_DIR")
  JOBID_MDR3=$(sbatch --parsable \
    -J rep3 \
    $PARTITION \
    --export=ALL \
    --dependency=afterok:${JOBID_SETUP} \
    "$SLURM_DIR/production/rep3/rep3.slurm" \
    "$SIMULATION_DIR")
  echo "JOBID_MDR1=${JOBID_MDR1}"
  echo "JOBID_MDR2=${JOBID_MDR2}"
  echo "JOBID_MDR3=${JOBID_MDR3}"

else
  echo "Error: Unknown SYSTEM type '$SYSTEM'" >&2
  exit 1
fi



