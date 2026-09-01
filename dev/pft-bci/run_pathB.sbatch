#!/bin/bash
#SBATCH --job-name=bci-pft-softmem
#SBATCH --output=bci-pft-softmem_%A_%a.out
#SBATCH --error=bci-pft-softmem_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpu_4,gpu_8,gpu_4_a100,gpu_4_h100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --array=1-3
################################################################################
# Path B: three species-resolution BCI conditions, one per array task.
#   1 = free (per-species baseline)   2 = ruger (fixed PFTs)   3 = learned (K=5)
#
# Prereqs on the cluster (see README.md):
#   - the FINN repo on branch soft-membership-pft (its FINN must be the version
#     installed in the container: finn_membership relies on FINN:::finn_class).
#   - dev/pft-bci/data/ present (run 01_build_species_obs.R once, or copy it).
#   - container ~/finn-r-torch.sif with FINN (this branch) + torch + CUDA.
#
# Submit from the FINN repo root:  sbatch dev/pft-bci/run_pathB.sbatch
################################################################################
export LANG=C.UTF-8; export LC_ALL=C.UTF-8; export R_PROFILE_USER=/dev/null

CONDS=(_ free ruger learned)                 # 1-indexed for SLURM_ARRAY_TASK_ID
export COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
export EPOCHS=8000
export PATCHES=25
export DEVICE=cuda
export BATCH=15
export LR=0.01
export SEED=42

echo "Path B condition: $COND  ($(date))"
singularity exec --nv ~/finn-r-torch.sif Rscript dev/pft-bci/03_pathB_fit.R
echo "done: $COND  ($(date))"
