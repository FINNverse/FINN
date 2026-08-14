# =============================================================================
# FINN | data-raw/nwfva/make_nwfva_gy.R
#
# Build FINN's bundled growth-&-yield reference dataset from the NW-FVA yield
# tables (Albert et al. 2021), so a first FINN version ships with real German
# growth-and-yield information for the five main commercial species.
#
# SOURCE (not redistributed here; install once to regenerate):
#   et.nwfva 0.1.0  --  Nuske, Staupendahl, Albert (NW-FVA), GPL (>= 2)
#   Zenodo record:  https://doi.org/10.5281/zenodo.7207597
#   Underlying data: Albert M., Nagel J., Schmidt M., Nagel R.-V., Spellmann H.
#     (2021) Eine neue Generation von Ertragstafeln fuer Eiche, Buche, Fichte,
#     Douglasie und Kiefer. Zenodo. https://doi.org/10.5281/zenodo.6343906
#   Tables: "<Art> gestaffelte Durchforstung 2021" (graduated thinning from
#   above is the underlying management concept).
#
# To regenerate:
#   install.packages("<path>/et.nwfva_0.1.0.tar.gz", repos = NULL, type = "source")
#   Rscript data-raw/nwfva/make_nwfva_gy.R
# =============================================================================

stopifnot(requireNamespace("et.nwfva", quietly = TRUE))
library(et.nwfva)

species     <- c("Eiche", "Buche", "Fichte", "Douglasie", "Kiefer")
yield_class <- -1:4   # relative Ertragsklasse; the tables interpolate/extrapolate

rows <- list()
for (s in species) for (b in yield_class) {
  d <- tryCatch(suppressWarnings(suppressMessages(et_tafel(s, bon = b))),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) next
  d$Art <- s
  rows[[paste(s, b)]] <- d
}
gy <- do.call(rbind, rows)
gy <- gy[order(gy$Art, gy$Ekl, gy$Alter), ]

# --- derived quantities FINN calibrates against -----------------------------
# Diameter increment of the remaining stand (cm/yr), per species x yield class.
gy$dDg_dt <- ave(gy$Dg,    gy$Art, gy$Ekl, FUN = function(x) c(NA, diff(x))) /
             ave(gy$Alter, gy$Art, gy$Ekl, FUN = function(x) c(NA, diff(x)))
# Thinning as the management signal: what fraction of stems / basal area is
# removed at each entry, and at what relative size (Dg_aus/Dg: <1 from below,
# >1 from above). These calibrate the management operator's thinning_intensity
# and thinning_size_bias (see dev/gy_growth_parameterisation.md).
gy$thin_N_frac <- gy$N_aus / (gy$N + gy$N_aus)
gy$thin_G_frac <- gy$G_aus / (gy$G + gy$G_aus)
gy$Dg_aus_rel  <- gy$Dg_aus / gy$Dg

# Column glossary (NW-FVA):
#   Ekl age N Hg H100 G Dg Dw V  = remaining stand after thinning
#     N stems/ha, Hg mean height m, H100 top height m, G basal area m2/ha,
#     Dg quadratic mean dbh cm, Dw dbh of the tree of mean basal area, V volume m3/ha
#   N_aus G_aus Dg_aus V_aus       = removed (thinning) stand this period
#   iV increment, GWL total yield performance, dGZ mean annual increment

utils::write.csv(gy, "data-raw/nwfva/nwfva_gy_tables.csv", row.names = FALSE)

# When promoted to a shipped dataset:
#   nwfva_gy <- gy; usethis::use_data(nwfva_gy, overwrite = TRUE)  # -> data/nwfva_gy.rda
message(sprintf("wrote %d rows, %d species, Ekl %s, age %d-%d",
                nrow(gy), length(unique(gy$Art)),
                paste(range(gy$Ekl), collapse = ".."),
                min(gy$Alter), max(gy$Alter)))
