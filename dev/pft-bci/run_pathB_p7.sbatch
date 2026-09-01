#!/bin/bash -l
#SBATCH --job-name=bci-pft-p7
#SBATCH --output=bci-pft-p7_%A_%a.out
#SBATCH --error=bci-pft-p7_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de
#SBATCH --mem=50G
#SBATCH --array=1-4
################################################################################
# Path B, PERIOD7 (5-yr step -- the published cadence), refit with FINN 0.2.0 as
# a reference for the annual-timestep runs. One FINN step per census interval,
# period_length=1, no regeneration_saturation (not needed at 5-yr step).
#   1 = free   2 = ruger   3 = learned   4 = pft5 (the paper's actual variant)
# ~5x fewer timesteps than annual, so short walltime.
################################################################################
export LANG=C.UTF-8; export LC_ALL=C.UTF-8
# do NOT set R_PROFILE_USER / R_LIBS_USER (see run_pathB.sbatch note)
module load compiler/gnu/14.2 numlib/mkl/2025.0 math/R/4.4.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
cd "$HOME/bci-pft" || exit 1
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH"; exit 1; }

CONDS=(_ free ruger learned pft5)
export COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
export PERIOD=7
export SAT=0
export TAG=_p7
export EPOCHS=4000
export PATCHES=25
export DEVICE=cpu
export BATCH=15
export LR=0.01
export SEED=42

echo "Path B P7: COND=$COND PERIOD=7 SAT=0  ($(date))"
Rscript dev/pft-bci/03_pathB_fit.R || { echo "FATAL: Rscript exited $? for COND=$COND"; exit 1; }
echo "done: $COND  ($(date))"
