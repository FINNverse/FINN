#!/bin/bash
#SBATCH --job-name=bci-pft-softmem
#SBATCH --output=bci-pft-softmem_%A_%a.out
#SBATCH --error=bci-pft-softmem_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpu_h100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yannek.kaeber@biom.uni-freiburg.de
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --array=1-6
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

# --- cluster R/torch environment (uc3: FINN 0.2.0 in ~/Rlib-finn020) ---------
module load math/R/4.4.2-mkl-2025.0-gnu-14.2
export R_LIBS_USER=$HOME/Rlib-finn020        # scripts also prepend this in-R
cd "$HOME/bci-pft" || exit 1                 # run root containing dev/pft-bci/
# -----------------------------------------------------------------------------

# conditions (1-indexed). Tasks 1-4: annual timestep WITH regeneration_saturation
# (the correct 0.2.0 config) across the four groupings. Tasks 5-6: saturation OFF,
# to isolate the effect of the annual-timestep regeneration fix.
CONDS=(_ free ruger learned pft5 learned pft5)
SATS=( _ 1    1     1       1    0       0)
export COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
export SAT=${SATS[$SLURM_ARRAY_TASK_ID]}
export TAG=$([ "$SAT" = "0" ] && echo "_nosat" || echo "")
export EPOCHS=8000
export PATCHES=25
export DEVICE=gpu
export BATCH=15
export LR=0.01
export SEED=42

echo "Path B: COND=$COND SAT=$SAT TAG=$TAG  ($(date))"
Rscript dev/pft-bci/03_pathB_fit.R
echo "done: $COND SAT=$SAT  ($(date))"
