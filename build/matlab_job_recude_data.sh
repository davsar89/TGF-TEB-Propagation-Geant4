#!/bin/bash
##SBATCH --account=NN9526K --qos=preproc
##SBATCH --account=NN9526K --qos=devel
#SBATCH --account=NN9526K --partition=bigmem
#SBATCH --job-name=postproc
#SBATCH --time=0-1:0:0
#SBATCH --nodes=1 --ntasks-per-node=32
#SBATCH --mem-per-cpu=3G

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

module restore system
module load MATLAB/2017a

cd $SLURM_SUBMIT_DIR
matlab -nodisplay -r "convert_to_mat_file_simple"
