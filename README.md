# AutoMD-SAXS
Automated molecular dynamics in GROMACS with post-processing and optional SAXS-based analysis.

#######include pic of pipeline##########


AutoMD-SAXS is an automated workflow for the setup, simulation and analysis of proteins, disorderd proteins, and protein-ligand systems. Built for use on HPC using the Slurm job scheduler.

## Installation
### Software requirements:
  - GROMACS (local install or module loaded)
  - Slurm job scheduler
  - Conda (all relevant packages are found in automdsaxs.yml). 
  - ATSAS 
      - Download only required for SAXS-based analysis 
      - Current analysis is optimised for ATSAS 3.0.4 (recommended)
   
### Setting up the Python environment

```bash
# Create a new conda environment
conda env create -f automdsaxs.ym
# Activate environment
conda activate automdsaxs
```


## Usage:

### Run files
  - Protein (.pdb)
  - SAXS data (.dat)

Note: Both files must be present within the AutoMD-SAXS directory

### Simulation setup
```bash
cd /path/to/AutoMD-SAXS 
```

```bash
sh simulation_setup.sh -p *Protein*.pdb --saxs *SAXS*.dat 
```
Outputs the directory ```bash *Protein*_simulation ```

### R

h
  - run_MD.sh 


    
## User workflow
  1. sh simulation_setup.sh -p \<protein\>.pdb -s \<saxs.dat\>
     - Flag -s is optional
     - Ensure .pdb and .dat file are present within the AutoMD-SAXS directory
     - Output is dir: '\<protein\>_simulation'
     - Creates '\<protein\>_simulation/configurations.txt' containing directory and parameter variables
    
  2. sh run_MD.sh <protein>_simulation
     - Calls dir 'slurms' to run MD jobs via Slurm

## Citation

If you use AutoMD-SAXS in your research...

```bibtex
@article{carbonara2025,
  title={},
  author={},
  journal={},
  year={2025},
  doi={},
  url={}
}

```
Shield:
This work is licensed under...
