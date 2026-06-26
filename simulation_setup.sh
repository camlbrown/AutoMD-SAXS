#!/bin/bash

echo "                _        __  __ _____         _____         __   __ _____ 
     /\        | |      |  \/  |  __ \       / ____|  /\    \ \ / // ____|
    /  \  _   _| |_ ___ | \  / | |  | |_____| (___   /  \    \ V /| (___  
   / /\ \| | | | __/ _ \| |\/| | |  | |______\___ \ / /\ \    > <  \___ \ 
  / ____ \ |_| | || (_) | |  | | |__| |      ____) / ____ \  / . \ ____) |
 /_/    \_\__,_|\__\___/|_|  |_|_____/      |_____/_/    \_\/_/ \_\_____/ 
"

echo "
     :.                                     =*************************=  
     : .                                    *=  +    -    -   + -    =*
     :  .                                   * =    .  +      ..-    = *
     :   .                                  *  =    . .  -  . -    = -*
     :    .                                 *+  = = = = = = = = = =   * 
     :      . .                             *   = -         -  .  = + *
     :        . . .                         *   =     +   @       =   *      
     :           .  . .                     *-  =@@@    @@@@@   @@@=  * 
     :             .  .  .  .      .        *  @@ @@   @@   @@ @@ @@- * 
     :                .  .  .  . . .  .     *   =  @@@@@  +  @@   =   *   
     :                 .  . . . .  . . .    * + =   @@@         - =  +*    
     :                    .   . .. . . .    *   = = = = = = = = = =   *   
     :                       .   . .   .    *  =    -.   -   ..-   =  *      
     :                             .   . .  * =. +    -.  -  .   .  = *
     :                                      *=      .  +       +     =*
     :::::::::::::::::::::::::::::::::::::  =*************************=
"



set -euo pipefail

# Help 
if [[ "$1" == "-help" || "$1" == "-h" ]]; then
  echo "Simulation Setup"
  echo ""
  echo "Prior to running simulation_setup.sh:"
  echo "  a) Make sure your protein (.pdb) is located within the AutoMD-SAXS directory."
  echo "  b) Make sure your SAXS data (.dat) is located within the AutoMD-SAXS directory."
  echo ""
  echo "Usage: ./simulation_setup.sh -p {protein.pdb} -s {saxs.dat}"
  echo "-s flag is optional"
  exit 0
fi

# Determine the directory this script lives in
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

PROTEIN_FILE=""
SAXS_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p)
      PROTEIN_FILE="$2"
      shift 2
      ;;
    -s)
      SAXS_FILE="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 -p <protein.pdb> [-s <saxs.dat>]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use -h or --help for usage instructions."
      exit 1
      ;;
  esac
done

if [[ -z "$PROTEIN_FILE" ]]; then
  echo "Error: Protein file (-p) must be provided." >&2
  exit 1
fi
if [[ "$PROTEIN_FILE" != *.pdb ]]; then
  echo "Error: Protein file must have a .pdb extension." >&2
  exit 1
fi

# Handle optional SAXS 
if [[ -z "$SAXS_FILE" ]]; then
  SAXS_FILE="None"
else
  if [[ "$SAXS_FILE" != *.dat ]]; then
    echo "Error: SAXS file must have a .dat extension." >&2
    exit 1
  fi
fi

PROTEIN_NAME="${PROTEIN_FILE%.pdb}"
export BASE_DIR PROTEIN_FILE PROTEIN_NAME SAXS_FILE

# Initial user notes 
echo "Welcome to your MD simulation setup!"
echo ""
echo "Protein file = $PROTEIN_FILE"
echo "SAXS file = $SAXS_FILE"
echo ""
echo "Please select your system type:"
echo ""
echo "Note: if your protein exhibits high levels of flexibility or contains intrinsically disordered regions, please select '2) Intrinsically Disordered Protein'" 
echo ""

# Directory structure 
SIMULATION_DIR="$BASE_DIR/${PROTEIN_NAME}_simulation"
export SIMULATION_DIR

# Top-level directories
declare -A STATIC_DIRS=(
  [FF_CONVERT]="$BASE_DIR/ff_convert"
  [FORCE_FIELD_DIR]="$BASE_DIR/ff_files"
  [SLURM_DIR]="$BASE_DIR/slurms"
  [MDP_DIR]="$BASE_DIR/mdp_files"
  [GMXLIB]="$BASE_DIR/ff_files"
)
for var in "${!STATIC_DIRS[@]}"; do
  mkdir -p "${STATIC_DIRS[$var]}"
  export "$var"="${STATIC_DIRS[$var]}"
done

# Subdirectories under $SIMULATION_DIR
STAGE_DIRS=(
  production
  production/rep1
  production/rep2
  production/rep3
  minim1
  minim2
  nvt
  npt
  pdb2gmx
  solvate
  genion
  SAXS
  ligand_setup
)

for stage in "${STAGE_DIRS[@]}"; do
  full="$SIMULATION_DIR/$stage"
  mkdir -p "$full"

  if [[ $stage =~ ^production/rep([1-3])$ ]]; then
    num=${BASH_REMATCH[1]}
    var="REPEAT_DIR${num}"
  elif [[ "$stage" == "ligand_setup" ]]; then
    var="LIGAND_SETUP"
  else
    var="${stage^^}_DIR"
    var="${var//\//_}"
  fi

  export "$var"="$full"
done

# Processed and frame-extraction dirs
for i in 1 2 3; do
  proc_var="PROCESSED_R${i}"
  proc_dir="$SIMULATION_DIR/production/rep${i}/processed"
  mkdir -p "$proc_dir"
  export "$proc_var"="$proc_dir"

  ext_var="R${i}_EXTRACT_FRAMES"
  ext_dir="$proc_dir/extract_frames"
  mkdir -p "$ext_dir"
  export "$ext_var"="$ext_dir"
done

# Use a per-job copy of the .mdp templates. run_MD.sh parameterises the
# production .mdp at run time; the legacy script edited the git-tracked template
# in $BASE_DIR/mdp_files in place. Copying here keeps the templates pristine, and
# because every Slurm script reads $MDP_DIR they transparently use this copy.
MDP_SOURCE="$BASE_DIR/mdp_files"
MDP_DIR="$SIMULATION_DIR/mdp_files"
mkdir -p "$MDP_DIR"
cp "$MDP_SOURCE"/*.mdp "$MDP_DIR"/
export MDP_DIR

# Copy SAXS data if provided
if [[ "$SAXS_FILE" != "None" ]]; then
  mkdir -p "$SAXS_DIR"
  cp "$SAXS_FILE" "$SAXS_DIR/" || {
    echo "Running without SAXS data." >&2
    exit 1
  }
fi

# MD SETUP

# GROMACS

config_file="$SIMULATION_DIR/configurations.txt"
read -p "Do you load GROMACS via 'module load'? (y/N) " load_mod
if [[ "$load_mod" =~ ^[Yy]$ ]]; then
  read -p "Which module (e.g. gromacs/2021.2/intel)? " gmx_load

  cat >> "$config_file" <<EOF
# GROMACS Module Management
GMX_MODULE="$gmx_load"
load_gmx()   { module load "\$GMX_MODULE"; }
unload_gmx() { module unload "\$GMX_MODULE"; }

EOF

else
  cat >> "$config_file" <<'EOF'
# No module loading  use local installation
GMX_MODULE=""
load_gmx()   { :; }
unload_gmx() { :; }

EOF

fi

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Python env location

read -p "Full path to the Python interpreter for the automdsaxs env (e.g. '/home/you/.conda/envs/automdsaxs/bin/python'): " python_cmd
python_cmd=${python_cmd:-python}

cat >> "$config_file" <<EOF
PYTHON_CMD="$python_cmd"
EOF

echo ""
echo "---------------------------------------------------------------------"
echo ""

# System setup

# 1) System type & (conditionally) force‐field choice
echo "Select system type:"
echo "  1) Protein"
echo "  2) Protein-ligand"
read -rp "Enter 1 or 2: " sys_choice

case "$sys_choice" in
  1)
    SYSTEM="Protein"
    export SYSTEM
    echo "You chose: $SYSTEM"
    echo ""

    # Now ask the force‐field for Protein
    echo "Select your force field:"
    echo "  1) amber14sb"
    echo "  2) charmm36m"
    read -rp "Enter 1 or 2: " ff_choice

    case "$ff_choice" in
      1)
        export FORCE_FIELD="$FORCE_FIELD_DIR/amber14sb"
        FF_CONVERT_SUBDIR="amber"
        ;;
      2)
        export FORCE_FIELD="$FORCE_FIELD_DIR/charmm36m"
        FF_CONVERT_SUBDIR="charmm"
        ;;
      *)
        echo "Invalid force‐field choice; please run again." >&2
        exit 1
        ;;
    esac
    ;;
  2)
    SYSTEM="Protein-ligand"
    export SYSTEM
    echo "You chose: $SYSTEM"
    echo ""

    # For Protein-ligand always use amber14sb
    export FORCE_FIELD="$FORCE_FIELD_DIR/amber14sb"
    FF_CONVERT_SUBDIR="amber"
    ;;
  *)
    echo "Invalid system type; please run again and pick 1 or 2." >&2
    exit 1
    ;;
esac

echo "Using force field: $FORCE_FIELD"
echo ""









#echo "Select system type:"
#echo "  1) Protein"
#echo "  2) Protein-ligand"
#read -rp "Enter 1 or 2: " sys_choice
#
#case "$sys_choice" in
#  1) SYSTEM="Protein"        ;;
#  2) SYSTEM="Protein-ligand" ;;
#  *)
#    echo "Invalid choice; please run again and pick 1 or 2." >&2
#    exit 1
#    ;;
#esac
#export SYSTEM
#echo "System type set to: $SYSTEM"
#echo ""


# Force-field choice 
#echo "Select your force field:"
#echo "  1) amber14sb"
#echo "  2) charmm36m"
#read -rp "Enter 1 or 2: " ff_choice

#case "$ff_choice" in
#  1)
#    export FORCE_FIELD="$FORCE_FIELD_DIR/amber14sb"
#    FF_CONVERT_SUBDIR="amber"
#    ;;
#  2)
#    export FORCE_FIELD="$FORCE_FIELD_DIR/charmm36m"
#    FF_CONVERT_SUBDIR="charmm"
#    ;;
#  *)
#    echo "Invalid choice; please run again and pick 1 or 2." >&2
#    exit 1
#    ;;
#esac

#echo "Using force field: $FORCE_FIELD"
#echo ""

# Protein Preparation Wizard check 
read -rp "Did you use the Protein Preparation Wizard to prepare your system for MD? (y/N) " ppw_ans
if [[ "$ppw_ans" =~ ^[Yy] ]]; then
  export PPW=yes
else
  export PPW=no
fi
echo "PPW set to: $PPW"

# Pick the right conversion directory root 
if [[ "$PPW" = "yes" ]]; then
  CONVERT_ROOT="$FF_CONVERT"
else
  CONVERT_ROOT="$FF_CONVERT/non_ppw_convert"
fi
echo "Will convert from: $CONVERT_ROOT/$FF_CONVERT_SUBDIR"
echo ""

# Run the appropriate conversion 
echo "Converting PDB with $FF_CONVERT_SUBDIR script in $CONVERT_ROOT..."
cp "$PROTEIN_FILE" "$CONVERT_ROOT/$FF_CONVERT_SUBDIR/"
pushd "$CONVERT_ROOT/$FF_CONVERT_SUBDIR" >/dev/null

if [[ "$FF_CONVERT_SUBDIR" = "amber" ]]; then
  sh amber_convert.sh "$PROTEIN_FILE"
else
  sh charmm_convert.sh "$PROTEIN_FILE"
fi

# Copy the generated GMX.pdb into your workflow
cp GMX.pdb "$SIMULATION_DIR"
cp GMX.pdb "$PDB2GMX_DIR"
cp GMX.pdb "$LIGAND_SETUP"
rm GMX.pdb


popd >/dev/null
echo "Conversion complete."
echo ""


# Box-shape 
echo "Now choose your box shape:"
echo "  1) Globular (dodecahedron)"
echo "  2) Anisotropic (rectangular)"
read -rp "Enter 1 or 2: " box_choice

case "$box_choice" in
  1) export BOX_SHAPE="dodecahedron" ;;
  2) export BOX_SHAPE="rectangular"  ;;
  *) echo "Invalid choice; exiting." >&2; exit 1 ;;
esac

echo "Box shape set to: $BOX_SHAPE"


echo ""
echo "---------------------------------------------------------------------"
echo ""

# Ionic strength

echo "Choose your ionic concentration (M) e.g. 0.15 is physiological salt concentration"
echo ""
echo "Note: if integrating SAXS data into the simulation, please use the experimental concentration"
echo ""

valid_input=false

while [ "$valid_input" = false ]; do
    read -p "Enter a number: " ionic_concentration

    if [[ $ionic_concentration =~ ^[0-9]*(\.[0-9]+)?$ && $(echo "$ionic_concentration > 0" | bc) -eq 1 ]]; then
        valid_input=true
    else
        echo "Error: Please enter a valid number greater than zero."
    fi
done

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Simulation time 

echo "Select simulation length."
echo ""

valid_input=false

while [ "$valid_input" = false ]; do
    read -p "Enter a number (nanoseconds): " simulation_time

    if [[ $simulation_time =~ ^[0-9]+$ && $simulation_time -gt 0 ]]; then
        valid_input=true
    else
        echo "Error: Please enter a valid non-zero number."
    fi
done

# Dmax

if [[ "$SAXS_FILE" != "None" ]]; then

  # Dmax
  echo "What is the maximum scattering dimension (Dmax) described by the experimental SAXS data?"
  echo ""
  echo "Note: if unknown please enter the largest dimension described by your protein model"
  echo ""

  valid_input=false

  while [ "$valid_input" = false ]; do
      read -p "Enter a number (nanometers): " dmax

      if [[ $dmax =~ ^[0-9]*(\.[0-9]+)?$ ]] && (( $(echo "$dmax > 0" | bc -l) )); then
          valid_input=true
      else
          echo "Error: Please enter a valid number greater than zero."
      fi
  done

else
  dmax="Model"
fi

if [[ "$SAXS_FILE" != "None" ]]; then
    cp $SAXS_FILE $SAXS_DIR
fi

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Disulfides 

echo "Does your system contain disulfide bonds?"
echo ""
valid_input=false

while [ "$valid_input" = false ]; do
    read -p "Enter 'yes' (y) or 'no' (n): " disulfide

    disulfide=$(echo "$disulfide" | tr '[:upper:]' '[:lower:]')

    if [[ "$disulfide" == "yes" || "$disulfide" == "y" || "$disulfide" == "no" || "$disulfide" == "n" ]]; then
        valid_input=true
    else
        echo "Error: Please enter either 'yes' (y) or 'no' (n)."
    fi
done


echo ""
echo "---------------------------------------------------------------------"
echo ""

# Mail

read -p "Would you like to be emailed upon simulation completion? (y/N) " mail_yn
if [[ "$mail_yn" =~ ^[Yy] ]]; then
  read -p "Enter your email address: " email_addr
else
  email_addr=""
fi

cat >> "$config_file" <<EOF
EMAIL_ADDR="$email_addr"
EOF

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Partition

read -p "Would you like to specify the SLURM partition now? (y/N) " partition_yn
if [[ "$partition_yn" =~ ^[Yy] ]]; then
  read -p "Enter the name of the partition you wish to run on: " partition_name
  PARTITION="--partition=${partition_name}"
else
  PARTITION=""
fi

cat >> "$config_file" <<EOF
PARTITION="${PARTITION}"
EOF

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Export variables to config file

touch $SIMULATION_DIR/monitor.txt
echo " " >> "$config_file"
echo "#Simulation variables " >> "$config_file"
echo "BASE_DIR=$BASE_DIR" >> "$config_file"
echo "FF_CONVERT=$FF_CONVERT" >> "$config_file"
echo "SIMULATION_DIR=$SIMULATION_DIR" >> "$config_file"
echo "PRODUCTION_DIR=$PRODUCTION_DIR" >> "$config_file"
echo "REPEAT_DIR1=$REPEAT_DIR1" >> "$config_file"
echo "REPEAT_DIR2=$REPEAT_DIR2" >> "$config_file"
echo "REPEAT_DIR3=$REPEAT_DIR3" >> "$config_file"
echo "PROCESSED_R1=$PROCESSED_R1" >> "$config_file"
echo "PROCESSED_R2=$PROCESSED_R2" >> "$config_file"
echo "PROCESSED_R3=$PROCESSED_R3" >> "$config_file"
echo "R1_EXTRACT_FRAMES=$R1_EXTRACT_FRAMES" >> "$config_file"
echo "R2_EXTRACT_FRAMES=$R2_EXTRACT_FRAMES" >> "$config_file"
echo "R3_EXTRACT_FRAMES=$R3_EXTRACT_FRAMES" >> "$config_file"
echo "MINIM1_DIR=$MINIM1_DIR" >> "$config_file"
echo "MINIM2_DIR=$MINIM2_DIR" >> "$config_file"
echo "NVT_DIR=$NVT_DIR" >> "$config_file"
echo "NPT_DIR=$NPT_DIR" >> "$config_file"
echo "PDB2GMX_DIR=$PDB2GMX_DIR" >> "$config_file"
echo "SOLVATE_DIR=$SOLVATE_DIR" >> "$config_file"
echo "GENION_DIR=$GENION_DIR" >> "$config_file"
echo "SAXS_DIR=$SAXS_DIR" >> "$config_file"
echo "FORCE_FIELD_DIR=$FORCE_FIELD_DIR" >> "$config_file"
echo "FORCE_FIELD=$FORCE_FIELD" >> "$config_file"
echo "SLURM_DIR=$SLURM_DIR" >> "$config_file"
echo "MDP_DIR=$MDP_DIR" >> "$config_file"
echo "BOX_SHAPE=$BOX_SHAPE" >> "$config_file"
echo "GMXLIB=$BASE_DIR/ff_files" >> "$config_file"
echo "PROTEIN_FILE=$PROTEIN_FILE" >> "$config_file"
#echo "SYSTEM=\"${options[$((choice-1))]}\"" >> "$config_file"
echo "SYSTEM=$SYSTEM" >> "$config_file"
echo "LIGAND_SETUP=$LIGAND_SETUP" >> "$config_file"
echo "IONIC_CONCENTRATION=$ionic_concentration" >> $config_file
echo "SIMULATION_TIME=$simulation_time" >> "$config_file"
echo "DISULFIDE=$disulfide" >> "$config_file"
echo "SAXS_FILE=$SAXS_FILE" >> "$config_file"
echo "DMAX=$dmax" >> "$config_file"
echo ""

# Validate the generated configuration with the Python package before the user
# runs run_MD.sh. Non-fatal: a preview only, so setup still succeeds if the
# package is unavailable.
echo "Validating configuration with automd_saxs..."
PDB_ARG=()
[ -f "$PDB2GMX_DIR/GMX.pdb" ] && PDB_ARG=(--pdb "$PDB2GMX_DIR/GMX.pdb")
PYTHONPATH="$BASE_DIR${PYTHONPATH:+:$PYTHONPATH}" "${python_cmd:-python3}" \
  -m automd_saxs plan \
  --config "$config_file" \
  --work-dir "$BASE_DIR" \
  --slurm-dir "$SLURM_DIR" \
  --mdp-dir "$MDP_DIR" \
  "${PDB_ARG[@]}" \
  || echo "(validation preview unavailable; ensure the automd_saxs package is importable)"
echo ""

cat <<EOM
Configurations file has been created: $config_file

Please check your system settings are correct in $SIMULATION_DIR/configurations.txt prior to executing 'run_MD.sh'

If you are happy with your configurations, run MD using:
"sh run_MD.sh <input_pdb>_simulation"  
EOM

