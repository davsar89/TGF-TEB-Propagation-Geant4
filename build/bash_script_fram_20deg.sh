#!/bin/bash
## Project:
#SBATCH --account=NN9526K -p normal
## Job name:
#SBATCH --job-name=TGF-TEB-propa
## Wall time limit:
#SBATCH --time=0-24:0:0
## Number of nodes and task er node
#SBATCH --nodes=8 --ntasks-per-node=32
#SBATCH --mail-user=david.sarria@uib.no

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default

module load foss/2017b
module load GCC/6.4.0-2.28
module load GCCcore/6.4.0
module load OpenMPI/2.1.1-GCC-6.4.0-2.28
module load gompi/2017b
module load Python/3.6.2-foss-2017b
module load Boost/1.66.0-foss-2018a-Python-2.7.14

module list             # List loaded modules, for easier debugging

## Run the application

cd $SLURM_SUBMIT_DIR

srun python run_on_multiple_cpu_20deg.py

wait

exit 0

