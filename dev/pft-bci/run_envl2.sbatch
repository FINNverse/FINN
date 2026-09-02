#!/bin/bash -l
#SBATCH --job-name=bci-envl2
#SBATCH --output=bci-envl2_%A_%a.out
#SBATCH --error=bci-envl2_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu
#SBATCH --mem=50G
#SBATCH --array=1-6
# env_l2 sweep: pft5 & learned x lambda {0, 0.05, 0.2}. Uses patched FINN (~/Rlib-finnfix).
module load compiler/gnu/14.2 numlib/mkl/2025.0 math/R/4.4.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
cd "$HOME/bci-pft" || exit 1
CONDS=(_ pft5 pft5 pft5 learned learned learned)
L2S=(  _ 0    0.05 0.2  0       0.05    0.2)
export COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
export ENV_L2=${L2S[$SLURM_ARRAY_TASK_ID]}
export EPOCHS=2000
export PATCHES=25
echo "envl2: COND=$COND ENV_L2=$ENV_L2 ($(date))"
Rscript dev/pft-bci/test_envl2.R || { echo "FATAL rc=$?"; exit 1; }
