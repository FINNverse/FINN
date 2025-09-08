# End-to-end: mainland USA polygon (no AK/HI/islands) + sample 1000 site points,
# using terra + geodata. CRS is taken from `temp` via crs(temp).
# Also creates a data.table `sites_dt` with a siteID column.

library(terra)
library(geodata)
library(data.table)

# 0) Target CRS from your existing object
tgt_crs <- crs(temp)

# 1) Countries from geodata (Natural Earth-like)
countries <- geodata::world(path = tempdir())     # WGS84

# 2) Subset to USA
usa <- countries[countries$GID_0 == "USA", ]

# 3) Clip to contiguous US via a bounding rectangle; drop Alaska/Hawaii/territories
conus_ext  <- ext(-125, -66.5, 24.5, 49.5)       # lon/lat bounds
conus_rect <- as.polygons(conus_ext, crs = crs(usa))

usa_mainland_ll <- aggregate(intersect(usa, conus_rect))  # dissolve to single polygon

# 4) Reproject mainland polygon to match `temp`/target CRS
usa_mainland <- project(usa_mainland_ll, tgt_crs)

# 5) SAMPLE 1000 site points inside mainland polygon (random sample)
#    (If you prefer regular lattice, use method="regular" instead.)
set.seed(42)  # reproducible
sites <- spatSample(usa_mainland, size = 1000, method = "regular")
sites
# 6) Keep only points strictly inside polygon (safety check)
sites <- sites[usa_mainland, ]

# 7) Make a data.table with a unique siteID
sites_df <- as.data.frame(sites, geom = "XY")  # WKT geometry column
sites_dt <- data.table(sites_df)
sites_dt[,siteID := 1:.N,]

# 8) Quick plot
plot(usa_mainland, col = "lightgray", border = "black")
points(sites, pch = 20, cex = 0.6, col = "red")

# Example: WorldClim BIO1 (temp) & BIO12 (prec)
temp <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
tempmax <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_5.tif")
tempmin <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_6.tif")
prec <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")
precseas <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_15.tif")
precwarmq <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_18.tif")

v_proj <- project(sites, crs(temp))
# env_dt[, temp := extract(temp, v_proj, method = "near")[, 2] ]
# env_dt[, prec := extract(prec, v_proj, method = "near")[, 2] ]
sites_dt[, ":="(
  temp = extract(temp, v_proj, method = "near")[, 2],
  tempmax = extract(tempmax, v_proj, method = "near")[, 2],
  tempmin = extract(tempmin, v_proj, method = "near")[, 2],
  prec = extract(prec, v_proj, method = "near")[, 2],
  precseas = extract(precseas, v_proj, method = "near")[, 2],
  precwarmq = extract(precwarmq, v_proj, method = "near")[, 2]
),]

envUS_dt <- sites_dt[,-c("x","y")]

# simulate ####
env_dt2 <- data.table()
for(i in 1:500){
  env_dt_temp <- envUS_dt
  env_dt_temp$year = i
  env_dt2 <- rbind(env_dt2, env_dt_temp)
}

for(i in names(scale_vals)){
  env_dt2[[i]] = (env_dt2[[i]]-scale_vals[[i]]$mean)/scale_vals[[i]]$sd
}

species_dt <- unique(inputs$obs_dt[,.(species, species_name)])
sim2 <- m1$simulate(env_dt2, init_cohort = NULL, device = "gpu", patches = 20, patch_size = 0.06)
p_dat <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]

wide_dt <- sim2$wide$site
wide_dt <- merge(wide_dt, species_dt, by = "species", all.x = TRUE)
wide_dt[,.(ba = mean(ba)), by = species_name][order(ba)]

library(ggplot2)
ggplot(p_dat[year <= 50],
       aes(x = year, y = value, color = factor(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")


