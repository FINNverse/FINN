#!/bin/bash -l
#SBATCH --job-name=bci-pft-validate
#SBATCH --output=bci-pft-validate_%j.out
#SBATCH --error=bci-pft-validate_%j.err
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu
#SBATCH --mem=50G
export LANG=C.UTF-8; export LC_ALL=C.UTF-8
module load compiler/gnu/14.2 numlib/mkl/2025.0 math/R/4.4.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
cd "$HOME/bci-pft" || exit 1
Rscript dev/pft-bci/05_validate.R
