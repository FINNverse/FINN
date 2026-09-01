#!/bin/bash
#SBATCH --job-name=bci-pft-softmem
#SBATCH --output=bci-pft-softmem_%A_%a.out
#SBATCH --error=bci-pft-softmem_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpu_4,gpu_8,gpu_4_a100,gpu_4_h100   # adjust to your uc3 allocation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --array=1-4
################################################################################
# Path B: four BCI conditions, one per array task.
#   1 = free (per-species baseline)      3 = learned (soft K=5)
#   2 = ruger (fixed 5 PFTs)             4 = pft5 (published 5-PFT data, this FINN)
#
# Runs FINN directly (no container). Prereqs:
#   - FINN installed in the cluster R env from branch soft-membership-pft
#     (finn_membership uses FINN:::finn_class). Current FINN is compatible with
#     the FINNetAl setup -- verified locally.
#   - dev/pft-bci/data/ present (run 01_build_species_obs.R once, or copy it).
#   - for pft5: the FINNetAl pft-period35-25patches CSVs at the path in 03.
#
# Submit from the FINN repo root:  sbatch dev/pft-bci/run_pathB.sbatch
################################################################################
export LANG=C.UTF-8; export LC_ALL=C.UTF-8; export R_PROFILE_USER=/dev/null

# --- your cluster R/torch environment (uncomment/adjust as needed) -----------
# module load devel/miniconda   ;  conda activate finn
# module load math/R/4.5.0
# -----------------------------------------------------------------------------

CONDS=(_ free ruger learned pft5)            # 1-indexed for SLURM_ARRAY_TASK_ID
export COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
export EPOCHS=8000
export PATCHES=25
export DEVICE=cuda
export BATCH=15
export LR=0.01
export SEED=42

echo "Path B condition: $COND  ($(date))"
Rscript dev/pft-bci/03_pathB_fit.R
echo "done: $COND  ($(date))"
