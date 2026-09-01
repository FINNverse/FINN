# FINN software note — working scaffold

`finn-software-note.Rmd` is a **structure-only scaffold** for the FINN software
note. Its section structure is taken from `FINN-true-software-note.md`
(Yannek's draft); its code chunks are the actual worked examples from the package
vignettes (`vignettes/C-Data_preparation.Rmd.orig`,
`vignettes/D-Fit_to_FIA.Rmd.orig`), reused directly. **No manuscript prose was
written** — every `<!-- PROSE ... -->` marks where Yannek's text belongs.

Style references for the writing: r3PG note (doi:10.1111/2041-210X.13474) and
cito note (doi:10.1111/ecog.07143).

## Precompile (like the vignettes)

The fits need a torch backend, so the note is precompiled once and the knitted
output committed (same pattern as `vignettes/build.R`). Chunks are `eval = FALSE`
by default so the scaffold knits for a structure check without training.

```r
# with a torch backend available:
# 1) set eval = TRUE in the setup chunk
# 2) knit once and commit the output + figures/
rmarkdown::render("paper/finn-software-note.Rmd")
```

## To-do to finish the manuscript

**Writing (Yannek — prose, not to be auto-generated):**
1. Fill each `<!-- PROSE -->` slot from `FINN-true-software-note.md` (Abstract,
   Introduction, workflow lead-in, the four case-study sections, Conclusion).
2. Confirm the author/affiliation/ORCID block and corresponding-author line.

**Analysis / build:**
3. Flip `eval = TRUE`, precompile with a torch backend, commit the knitted
   output and `figures/` (mirror `vignettes/build.R`).
4. Decide which figures to include: the concept figure (`Concept_Figure.png` /
   `vignettes/FINN-overview.jpg`), the growth-ALE comparison (mechanistic vs NN),
   and the predicted-vs-observed assessment. Add the concept figure file.
5. Verify every code chunk against the **CRAN 0.2.0** API (they are copied from
   the shipped 0.2.0 vignettes, so they should match as-is).

**Manuscript infrastructure:**
6. Choose the target journal + template (the framing and the two style refs point
   to MEE; consider `rticles` or the journal's Rmd/Quarto template) and switch the
   YAML `output:` accordingly.
7. Bibliography: turn `references.md` (from the draft) into a `.bib` + a journal
   CSL, and wire in `[@key]` citations.
8. Author contributions, data-availability (CRAN + Zenodo/FIA), and a cover
   letter / submission checklist for the target journal.
