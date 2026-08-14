# Parameterising FINN from NW-FVA yield tables + PROFOUND

*A first FINN version that ships with real growth-and-yield information and can
simulate management. This note records the data, the plan, and — explicitly —
what the data cannot tell us.*

## 1. What we now have

**NW-FVA yield tables** (`et.nwfva` 0.1.0, Albert et al. 2021; Zenodo
[7207597](https://doi.org/10.5281/zenodo.7207597)). Five main species — Eiche,
Buche, Fichte, Douglasie, Kiefer — as *"gestaffelte Durchforstung 2021"* tables.
Extracted to `data-raw/nwfva/nwfva_gy_tables.csv` (642 rows; yield classes -1..4;
age 15–180). Per species × yield class × 5-yr age step:

| block | columns | meaning |
|---|---|---|
| remaining stand | `N, Hg, H100, G, Dg, Dw, V` | after thinning |
| **removed (thinning)** | `N_aus, G_aus, Dg_aus, V_aus` | the management signal |
| performance | `iV, GWL, dGZ` | increment / total yield |
| derived here | `dDg_dt, thin_N_frac, thin_G_frac, Dg_aus_rel` | FINN calibration targets |

**Schätzhilfen** (`schaetzhilfen.pdf`): site-index lookup, Alter × Höhe → yield
class (dGZ100). The classifier that places a real stand on the right table row.

**PROFOUND** (Bílý Kříž): 19 annual censuses of an *unmanaged* Norway-spruce
stand (stems 602→313). Complementary — see §4.

## 2. Growth: yield tables → FINN's `growth`

The tables give the diameter trajectory `Dg(age)` per species and site index, i.e.
increment `dDg_dt` as a function of size and productivity. Example (Fichte, Ekl 1):
`dDg_dt` peaks ~0.74 cm/yr around Dg≈15 cm and declines to ~0.34 cm/yr at Dg≈45 cm
— the classic decelerating individual-diameter growth FINN's `growth` process is
meant to express.

**Plan.** Treat each (species × yield class) trajectory as a managed-stand growth
curve and calibrate FINN's `growth` so simulated cohorts reproduce `Dg(age)` across
the site-index fan. Site index (from the Schätzhilfen / `et_bonitaet`) enters as the
environmental gradient the growth process reads. This gives FINN real, regionally
valid growth for the five species out of the box.

## 3. Management: the Durchforstungshilfe **does** map onto the module

Question posed: can the FVA Durchforstungshilfe / Z-Baum logic be represented by the
parametric management operator? **Yes — and the yield tables supply the calibration
target.** The removal columns show the graduated regime directly (Fichte, Ekl 1):

| age | Dg | `thin_G_frac` | `Dg_aus_rel` | reading |
|---:|---:|---:|---:|---|
| 25 | 9.7 | 0.25 | 0.80 | heavy tending **from below** |
| 35 | 17.0 | 0.16 | 0.69 | release, from below |
| 50 | 24.4 | 0.10 | 1.24 | **from above** (harvest the strong) |
| 55 | 26.0 | 0.09 | 1.33 | from above |
| 90 | 38.6 | 0.06 | 0.97 | near-neutral maintenance |

So the NW-FVA management concept is: **strong early Niederdurchforstung → mid-life
Hochdurchforstung → light neutral maintenance**, with intensity falling
monotonically (0.25→0.05). Mapped onto the operator (see the WET2024 catalogue):

- `thinning_intensity` (Durchforstungsstärke) = `thin_G_frac`, size-declining.
- `thinning_size_bias` (Durchforstungsart) = a monotone function of `Dg_aus_rel`;
  it is **not constant** — it swings negative→positive with size. **Design
  refinement:** promote `thinning_size_bias` from a scalar to a function of relative
  size (or age), fitted to `Dg_aus_rel`. This is exactly the Durchforstungshilfe's
  "hold the Z-Bäume by removing their competitors, shifting up the canopy as the
  stand matures".
- Z-Baum selection itself (individual, spatial, marked crop trees) remains
  approximated by the size signature, not reproduced tree-by-tree — the standing
  limitation already stated in the catalogue.

Net: the Durchforstungshilfe/Z-Baum regime is representable, and calibratable to the
`*_aus` columns. The FVA-BW interactive Durchforstungshilfe (h/d targets, Z-Baum
counts per Oberhöhe) can refine the exact `thinning_intensity(size)` schedule if we
want the BW numbers rather than the NW-FVA ones.

## 4. The joint dataset — and handling missing processes **explicitly**

The plan is one combined calibration set. But the two sources inform *different*
FINN processes, and neither covers all four. Stating this plainly is the point:

| FINN process | NW-FVA yield tables | PROFOUND (Bílý Kříž) | decision |
|---|---|---|---|
| **growth** | ✔ strong, all 5 species × site index | ✔ one spruce stand | fit primarily to yield tables; PROFOUND as an independent check |
| **management/thinning** | ✔ explicit `*_aus` schedule | ✘ (unmanaged) | fit the operator to the yield-table removals |
| **mortality** (natural) | ✘ folded into thinning; no natural-death column | ✔ real self-thinning (602→313, unmanaged) | fit mortality to PROFOUND; yield tables give ~0 natural mortality in the managed, tended regime |
| **regeneration** | ✘ even-aged, starts at establishment N | ✘ planted, single cohort | **not identifiable from either** — see below |

**Explicit gap handling (do not fit silently):**

1. **Mortality.** Yield tables carry no natural mortality (removals are harvest, not
   death). Fitting growth to them alone would learn "no mortality". So mortality is
   fit to PROFOUND's unmanaged self-thinning and/or held at a documented prior;
   during yield-table fitting the mortality loss weight is set to ~0 (managed stands
   have their density set by thinning, not death). The combined fit must **weight
   mortality only where it is observed**, or it will be biased by construction — the
   same failure mode the BWI work hit (`dev`/FINN-bwi harvest-process note).
2. **Regeneration.** Neither source constrains recruitment: yield tables begin at a
   given establishment density, PROFOUND is a planted monoculture. Regeneration is
   therefore **fixed to a prior / switched off** for this first version and flagged
   as an assumption, not calibrated. `simulate_managed()` with `planting_rate`
   (Verjüngung/Vorbau) supplies establishment prescriptively instead of relying on a
   fitted regeneration process. This is stated wherever managed runs are reported.
3. **Off-support extrapolation.** Yield-table growth is *managed-stand* growth
   (tended density, thinning from above). Using it to predict unmanaged old-growth is
   off-support; the PROFOUND fit is the unmanaged anchor. Report which regime a run
   sits in.

## 5. Build order (proposed)

1. **Data** ✔ — `data-raw/nwfva/make_nwfva_gy.R` → `nwfva_gy_tables.csv`. Promote to
   a shipped `data/nwfva_gy.rda` once the schema is fixed (adds "GY tables ship with
   FINN"). `et.nwfva` stays a `Suggests` for regeneration; the derived table is
   vendored with attribution.
2. **Growth fit** — calibrate FINN `growth` to `Dg(age)` across species × site index;
   validate against PROFOUND spruce.
3. **Management** — make `thinning_size_bias` size-dependent; fit `thinning_intensity(size)`
   and the bias to the `*_aus` columns; check `simulate_managed()` reproduces a yield
   table's own thinning schedule.
4. **Mortality/regeneration** — mortality from PROFOUND (weighted where observed);
   regeneration fixed/prescriptive with the assumption documented.

## Decisions (locked 2026-08-14)

1. **Ship the tables in the package as an open, non-binary CSV** —
   `inst/extdata/nwfva_gy_tables.csv`, loaded via `nwfva_yield_tables()`. Not an
   `.rda`. Provenance/license beside it in `nwfva_gy_tables_SOURCE.md`.
2. **`thinning_size_bias` becomes a function of size, not a scalar** — required to
   capture the below→above swing (§3). v1 fits `bias(relative size)` to `Dg_aus_rel`.
3. **Site index enters as relative Ertragsklasse (Ekl) for v1.** *Flagged limitation:*
   Ekl is a productivity index, **not** a climatic response — using it as the growth
   "environment" conflates site quality with climate. A later version must replace/augment
   it with actual climate drivers and **critically assess the climatic response** of the
   fitted growth (this is a known v1 shortcut, not the end state).
4. **Use both thinning schedules** — the NW-FVA `*_aus` outcome *and* the FVA-BW
   Durchforstungshilfe decision rule (see §3).

## Licensing & citation (done)

- **Data:** Albert, M., Nagel, J., Schmidt, M., Nagel, R.-V., Spellmann, H. (2021)
  *Eine neue Generation von Ertragstafeln …* \[Dataset\], Zenodo,
  doi:10.5281/zenodo.6343906 (concept doi:10.5281/zenodo.14228056). **CC-BY-4.0**,
  © NW-FVA. Redistributable in a GPL package with attribution + "indicate changes";
  both satisfied in `nwfva_gy_tables_SOURCE.md` and the `nwfva_yield_tables()` docs.
  The data keeps its CC-BY-4.0 license, separate from FINN's GPL-3.
- **Extraction tool:** Nuske, Staupendahl, Albert (2022) *et.nwfva* 0.1.0, Zenodo
  doi:10.5281/zenodo.7207597, GPL(>=2). A `Suggests`, used only to regenerate the CSV.
- **Durchforstungshilfe:** Klädtke, J. (2024) *Durchforstungshilfe 2024*, FVA
  PRAXISNAH Heft 2. Used as design input (rule + numbers below), not redistributed.
- **TODO before any release:** add a data-source note to `DESCRIPTION` so the CC-BY-4.0
  provenance is visible in package metadata.

## Appendix: FVA-BW Durchforstungshilfe 2024 (Df-24) — the decision rule

Z-Baum-oriented Auslesedurchforstung, steered by **Oberhöhe** and the **h/d of the
Z-Bäume**:
- Below ~12 m Oberhöhe: select and mark Z-Bäume. 12–25 m: the release phase. Above
  25 m: thinning ends.
- **Stability target `h/d ≤ 75`** for Z-Bäume; if `h/d > 75`, remove at least one
  competitor (Bedränger) per Z-Baum.
- **Turnus is height-growth-driven:** re-decide after ~**3 m** Oberhöhen-Zuwachs (not a
  fixed number of years) — maps to `entry_interval_years` via the height curve.
- Hold Z-Bäume on a **Z-Baum-Norm** diameter trajectory by removing their strongest
  competitors → selective **from above** among non-crop trees; lagging diameter
  (red curves) ⇒ higher `thinning_intensity`.
- Max Z-Baum count per ha decreases with target diameter (Tab. 1; a graphic — exact
  NZB/Ziel-BHD numbers to be transcribed if we adopt the BW crop-tree densities).

This is the *rule* that produces the NW-FVA `*_aus` *outcome*: the two are
complementary — Df-24 for the prescriptive crop-tree logic, the yield tables for the
empirical stand-level removal to calibrate against.
