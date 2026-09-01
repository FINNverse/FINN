# Soft-membership PFTs on BCI (0.3.0 exploration)

Tests whether a learned soft membership lets FINN fit all BCI species at once
while borrowing statistical strength for rare species — and whether the learned
groups resemble Rüger et al. (2020)'s 5 PFTs. Nothing here touches the shipped
package; it all lives on branch `soft-membership-pft`.

## The idea
Each species' mechanistic parameters become `par_s = forward(A[s,] %*% Proto)`,
with **one shared membership** `A = softmax(logits)` `[N_species × K]` and
per-process prototypes `Proto` `[K × n_par]`. `finn_membership.R` implements this
as a subclass of `FINN:::finn_class`, overriding only the four `par_<type>`
active bindings — the whole FINN forward/`fit`/`predict` is reused unchanged.

## Files
| file | what |
|---|---|
| `finn_membership.R` | the subclass (learned or frozen membership) |
| `test_finn_membership.R` | end-to-end smoke test on the FIA example |
| `01_build_species_obs.R` | build species-level BCI obs/env/cohorts from FINNetAl |
| `02_pathA_membership.R` | Path A: membership on the regression surrogate (local) |
| `03_pathB_fit.R` | Path B: fit ONE condition of the full FINN experiment |
| `04_pathB_analyze.R` | Path B: rare-species gain + Rüger match |
| `run_pathB.sbatch` | bwUniCluster array (one condition per task) |

`data/` (git-ignored — derived from Dryad-licensed BCI) is produced by
`01_build_species_obs.R`, which reads the cleaned per-tree file from a sibling
`FINNetAl` checkout.

## Path B — the experiment
Four models, identical except for the grouping of the per-species demographic
parameters (the first three via `finn_membership` at species resolution):

- **free** — frozen identity membership → free per-species (baseline).
- **ruger** — frozen Rüger PFT one-hot → parameters tied by the 5 PFTs.
- **learned** — soft membership, `K = 5`, discovered by the fit.
- **pft5** — the *published* approach: data aggregated to 5 Rüger PFTs, refit
  with the current FINN as a same-version anchor.

A fixed 5-of-20 site hold-out measures prediction where training data is scarce.
**Primary** result: held-out error by rare vs common species (`04`).
**Secondary**: does `learned`'s membership match Rüger (ARI + confusion, `04`).

## Run it
Local (fast smoke): `COND=learned EPOCHS=20 PATCHES=4 Rscript dev/pft-bci/03_pathB_fit.R`

Cluster (real, ~overnight), from the FINN repo root:
1. `Rscript dev/pft-bci/01_build_species_obs.R`  (or copy `data/` over)
2. `sbatch dev/pft-bci/run_pathB.sbatch`   (array 1-4: free / ruger / learned / pft5)
3. back home: `Rscript dev/pft-bci/04_pathB_analyze.R`

**FINN version:** runs FINN directly (no container). The cluster R env must have
FINN installed from **this branch** (`finn_membership` uses `FINN:::finn_class`).
Current FINN is compatible with the FINNetAl setup — verified locally by refitting
both the 5-PFT and species data.

## Config notes
- `batchsize` is clamped to the number of training sites (FINN drops the last
  incomplete site-batch, so `batchsize > n_sites` → zero batches → no fit).
- Env effects are left per-species (free) in all three conditions, so the
  comparison isolates the pooling of the *demographic* parameters.
