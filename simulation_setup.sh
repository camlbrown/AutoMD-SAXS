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

# Validate protein input
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

# Copy SAXS data if provided
if [[ "$SAXS_FILE" != "None" ]]; then
  mkdir -p "$SAXS_DIR"
  cp "$SAXS_FILE" "$SAXS_DIR/" || {
    echo "Running without SAXS data." >&2
    exit 1
  }
fi

#---------------MD SETUP---------------

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

options=("Protein" "Intrinsically Disordered Protein" "Protein-ligand")

for i in "${!options[@]}"; do
  echo "$((i+1)). ${options[i]}"
done

echo ""

read -p "Enter the number corresponding to your system type: " choice

case $choice in
  1)
    echo "You chose Protein."
    echo "Amber14ffsb force field selected. If you would like to change force field see the manual"
    while true; do
      echo ""
      echo "--------------------------------------------------------------------------------------------"
      echo ""
      echo "Choose box shape:"
      echo ""
      echo "1. Globular"
      echo "2. Anisotropic"
      echo ""
      read -p "Enter the number corresponding to the box shape: " box_choice
      echo ""
      case $box_choice in
        1)
          export BOX_SHAPE="dodecahedron"
          break
          ;;
        2)
          export BOX_SHAPE="rectangular"
          break
          ;;
        *)
          echo "Invalid choice. Please choose again."
          ;;
      esac
    done
    export FORCE_FIELD="$FORCE_FIELD_DIR/amber14sb"
    cp $PROTEIN_FILE $FF_CONVERT/amber && cd $FF_CONVERT/amber
    sh amber_convert.sh $PROTEIN_FILE
    cp GMX.pdb $SIMULATION_DIR 
    cp GMX.pdb $PDB2GMX_DIR 
    rm GMX.pdb && cd $BASE_DIR  
    ;;
  2)
    echo "You chose Intrinsically disordered protein."
    echo "Amber14ffsb force field selected. If you would like to change force field see the manual"
    while true; do
      echo ""
      echo "--------------------------------------------------------------------------------------------"
      echo ""
      echo "Choose box shape:"
      echo ""
      echo "1. Globular"
      echo "2. Anisotropic"
      echo ""
      read -p "Enter the number corresponding to the box shape: " box_choice
      echo ""
      case $box_choice in
        1)
          export BOX_SHAPE="dodecahedron"
          break
          ;;
        2)
          export BOX_SHAPE="rectangular"
          break
          ;;
        *)
          echo "Invalid choice. Please choose again."
          ;;
      esac
    done
    cp $PROTEIN_FILE $FF_CONVERT/charmm && cd $FF_CONVERT/charmm
    sh charmm_convert.sh $PROTEIN_FILE
    cp GMX.pdb $SIMULATION_DIR
    cp GMX.pdb $PDB2GMX_DIR 
    rm GMX.pdb && cd $BASE_DIR  
    export FORCE_FIELD="$FORCE_FIELD_DIR/charmm36m"
    ;;

  3)
    echo "You chose Protein-ligand."
    while true; do
      echo "Choose box shape:"
      echo "1. Globular"
      echo "2. Anisotropic"
      read -p "Enter the number corresponding to the box shape: " box_choice
      case $box_choice in
        1)
          export BOX_SHAPE="octahedron"
          break
          ;;
        2)
          export BOX_SHAPE="triclinic"
          break
          ;;
        *)
          echo "Invalid choice. Please choose again."
          ;;
      esac
    done
    cp $PROTEIN_FILE $FF_CONVERT/amber && cd $FF_CONVERT/amber
    sh amber_convert.sh $PROTEIN_FILE
    cp GMX.pdb $SIMULATION_DIR 
    cp GMX.pdb $LIGAND_SETUP 
    cp GMX.pdb $PDB2GMX_DIR
    rm GMX.pdb && cd $BASE_DIR  
    export FORCE_FIELD="$FORCE_FIELD_DIR/amber14sb"
    ;;
  *)
echo "Invalid choice. Please select a valid number from the options provided."
    exit 1
    ;;
esac

echo ""
echo "---------------------------------------------------------------------"
echo ""

# Ionic strength

echo "Choose your ionic concentration (mM) e.g. 0.15 "
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
echo "SYSTEM=\"${options[$((choice-1))]}\"" >> "$config_file"
echo "LIGAND_SETUP=$LIGAND_SETUP" >> "$config_file"
echo "IONIC_CONCENTRATION=$ionic_concentration" >> $config_file
echo "SIMULATION_TIME=$simulation_time" >> "$config_file"
echo "DISULFIDE=$disulfide" >> "$config_file"
echo "SAXS_FILE=$SAXS_FILE" >> "$config_file"
echo "DMAX=$dmax" >> "$config_file"
echo ""

cat <<EOM
Configurations file has been created: $config_file

Please check your system settings are correct in $SIMULATION_DIR/configurations.txt prior to executing 'run_MD.sh'

If you are happy with your configurations, run MD using:
"sh run_MD.sh <input_pdb>_simulation"  
EOM

