# AutoMD-SAXS
Automated molecular dynamics in GROMACS with post-processing and optional SAXS-based analysis.

AutoMD-SAXS is an automated workflow for the setup, simulation and analysis of proteins, disorderd proteins, and protein-ligand systems. Integrated with the Slurm scheduler, this pipeline is designed for HPC.

**Software requirements:**
  - GROMACS 
  - Slurm job scheduler
  - Conda (all relevant packages are found in automdsaxs.yml)
  - ATSAS
      - Current analysis is optimised for ATSAS 3.0.4 (recommended)

**Scripts:**
  - simulation_setup.sh
  - run_MD.sh 

**Run requirements**
  - Protein file (.pdb)
  - SAXS data (.dat)

    
**User workflow**
  1. sh simulation_setup.sh -p <protein>.pdb -s <saxs.dat>
     - Flag -s is optional
     - Ensure .pdb and .dat file are present within the AutoMD-SAXS directory
     - output is dir: _<protein>_simulation_
     - creates configurations.txt containing directory and parameter variables
    
  2. sh run_MD.sh <protein>_simulation
     - calls dir _slurms_ to run MD jobs via Slurm

