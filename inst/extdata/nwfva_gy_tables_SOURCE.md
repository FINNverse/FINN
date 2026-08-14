# `nwfva_gy_tables.csv` — source, license, and attribution

This dataset is a **derived, reformatted extract** of the NW-FVA yield tables. It
is redistributed with FINN under the terms of its original license (below), kept
separate from FINN's own GPL-3 license.

## Original data (the source of these values)

> Albert, M., Nagel, J., Schmidt, M., Nagel, R.-V., Spellmann, H. (2021):
> *Eine neue Generation von Ertragstafeln für Eiche, Buche, Fichte, Douglasie und
> Kiefer* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.6343906
> (concept DOI, resolves to the latest version: https://doi.org/10.5281/zenodo.14228056)

**License: Creative Commons Attribution 4.0 International (CC-BY-4.0)** —
https://creativecommons.org/licenses/by/4.0/
Copyright: Nordwestdeutsche Forstliche Versuchsanstalt (NW-FVA).

CC-BY-4.0 permits redistribution and reuse (including in a GPL-licensed package)
provided the creators are credited and changes are indicated. Both are done here.

## Extraction tool

Values were read out of the R package that packages the above dataset:

> Nuske, R., Staupendahl, K., Albert, M. (2022): *et.nwfva: Forest Yield Tables for
> Northwest Germany and their Application.* R package version 0.1.0. Zenodo.
> https://doi.org/10.5281/zenodo.7207597 — License: GPL (>= 2).

## Changes made to the original (CC-BY "indicate changes" requirement)

- Extracted the five species tables ("… gestaffelte Durchforstung 2021") for
  relative yield classes (Ertragsklasse) −1 to 4 over the tabulated age range,
  via `et.nwfva::et_tafel()`, into one tidy long table (column `Art` added).
- Added four **derived** columns not in the original tables:
  `dDg_dt` (diameter increment, cm yr⁻¹, from consecutive `Dg`), `thin_N_frac`,
  `thin_G_frac` (removed stems / basal area as a fraction of the pre-thinning
  stand) and `Dg_aus_rel` (`Dg_aus / Dg`).
- No original tabulated value was altered; rows are the native 5-year age steps.

Regeneration script: `data-raw/nwfva/make_nwfva_gy.R`.

## Column glossary

`Art` species · `Ekl` relative yield class · `Alter` age (yr). Remaining stand
after thinning: `N` stems ha⁻¹, `Hg` mean height (m), `H100` top height (m),
`G` basal area (m² ha⁻¹), `Dg` quadratic mean dbh (cm), `Dw` dbh of the mean-basal-
area tree (cm), `V` volume (m³ ha⁻¹). Removed (thinning) stand: `N_aus, G_aus,
Dg_aus, V_aus`. Performance: `iV` increment, `GWL` total yield, `dGZ` mean annual
increment. Derived: `dDg_dt, thin_N_frac, thin_G_frac, Dg_aus_rel`.
