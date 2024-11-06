library(data.table)

out_dir = "data/calibration-data/BCI-1h-patch-V2"
if(!dir.exists(out_dir)) {
  dir.create(out_dir,recursive = T)
}

census_years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

total_area = 1000*500
area_of_square = total_area/500
sqrt(1000)
# BCI data was downloaded from https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
dt_list <- list()
for (i in 1:8) {
  # load rdata file
  load(paste0("data/calibration-data/BCI/bci.tree/bci.tree",i,".rdata"))
  dt <- data.table(get(paste0("bci.tree",i)))
  dt$census = census_years[i]
  dt_list[[i]] <- dt
}
load(paste0("data/calibration-data/BCI/bci.spptable.rdata"))
spptable <- data.table(bci.spptable)

# donwloaded from
# https://www.osti.gov/biblio/1771850
meteo_dt <- fread("data/calibration-data/BCI/BCI_1985_2018c_mod_2018substituted_20210129025453.csv")

cor(as.matrix(meteo_dt[,.(Temp_deg_C, RH_prc, SR_W_m2, Rainfall_mm_hr)]), method = "spearman")

# convert 12/31/84 20:00 to year, month, day, hour columns
meteo_dt[, c("date", "hour") := tstrsplit(DateTime, " ", fixed = TRUE)]
meteo_dt[, c("month", "day", "year") := lapply(tstrsplit(date, "/", fixed = TRUE), as.integer)]
meteo_dt[year > 80, year := as.integer(paste0(19,year)),]
meteo_dt[year < 80 & year >= 10, year := as.integer(paste0(20,year)),]
meteo_dt[year < 80 & year < 10, year := as.integer(paste0(200,year)),]
meteo_dt[, c("hour") := as.integer(gsub(":00", "", hour))]

env_dt = meteo_dt[,.(
  Prec = sum(Rainfall_mm_hr),
  SR_W_m2 = sum(SR_W_m2),
  Temp = mean(Temp_deg_C),
  RH_prc = mean(RH_prc)
), by = .(year, month, day)]

env_dt <- env_dt[,.(
  Prec = sum(Prec),
  SR_kW_m2 = sum(SR_W_m2)/1000,
  RH_prc = mean(RH_prc),
  T_mean = mean(Temp),
  T_max = max(Temp),
  T_min = min(Temp)
), by = .(year)]

# fwrite(env_dt, "data/calibration-data/BCI/env_dt.csv")
fwrite(env_dt[year >= 1985], paste0(out_dir, "/env_dt.csv"))
# # Load the GIFT library
# library(GIFT)
#
# # Get species data from the GIFT database
# species_data <- GIFT_species()
# species_data <- data.table(species_data)
#
# # Get metadata for available traits
# traits_meta <- data.table(GIFT_traits_meta())
#
# # View the metadata to identify the correct trait for "growth form"
# print(traits_meta)
#
# # Get trait data from the GIFT database
# traits_data <- GIFT_traits("1.2.1")
# traits_data <- data.table(traits_data)
#
# spp_dt <- merge(
#   spptable[,.(sp, species = paste(Genus,Species))],
#   traits_data[,.(trait_value_1.2.1, species = work_species, genus = )], by = "species")
#
# # Filter for tree species and select unique families
# tree_traits <- traits_data[traits_data$trait_value == "tree", ]
#
# # Extract unique family names from the filtered data
# tree_families <- unique(tree_traits$family)
#
# # Print the list of tree families
# print(tree_families)
#
# # Filter for tree species and select unique families
# tree_families <- unique(species_data[growth_form == "tree", "family"])
#
# # Print the list of tree families
# print(tree_families)

all_trees <- rbindlist(dt_list, fill=TRUE)
unique(all_trees[,.(treeID, census)])


all_trees <- merge(all_trees, spptable[,.(sp,Genus)], by = "sp")
uniqueN(all_trees$Genus)
all_trees <- all_trees[order(census)]
table(all_trees$status)

all_trees[status %in% c("A", "AD", "AR"), status := "A",]
all_trees[,species := as.integer(as.factor(Genus)),]
all_trees[, dbh_cm := dbh/10]
all_trees[, ":="(
  status_before = data.table::shift(status, 1, type = "lag"),
  dbh_cm_before = data.table::shift(dbh_cm, 1, type = "lag"),
  period_length = census - data.table::shift(census, 1, type = "lag")
  ),
  by = .(treeID)]


all_trees[, g := (dbh_cm-dbh_cm_before)/period_length,]

# g_max <- quantile(all_trees$g,0.9999, na.rm = T)
g_max = 5
all_trees[g > g_max | g < 0, g := NA]
# hist(all_trees$g, breaks = 100)


uniqueN(all_trees$quadrat)
x_length = 40
y_length = 25
# x_length = 100
# y_length = 100
all_trees <- all_trees[!is.na(gx) & !is.na(gy)]
all_trees[,":="(
  x_class = cut(gx, breaks = seq(0, 1000, x_length), labels = 1:length(seq(x_length, 1000, x_length)), include.lowest = T),
  y_class = cut(gy, breaks = seq(0, 500, y_length), labels = 1:length(seq(y_length, 500, y_length)), include.lowest = T)
),]
all_trees[is.na(x_class)]

all_trees[, patch := as.integer(as.factor(paste0(x_class,"-",y_class))),]
uniqueN(all_trees$treeID)
uniqueN(all_trees$patch)

# library(ggplot2)
# library(scales)
#
# # Generate random colors for the unique patches
# set.seed(123)  # Set seed for reproducibility
# n_patches <- length(unique(all_trees$patch))  # Get the number of unique patches
# random_colors <- hue_pal()(n_patches)  # Generate n random colors using hue palette

# ggplot(all_trees[census == 2015], aes(x = gx, y = gy, colour = factor(patch, ordered = FALSE))) +
#   geom_point(size = 0.01, alpha = 0.5) +
#   scale_color_manual(values = sample(random_colors)) +  # Assign random colors
#   theme_minimal() +
#   theme(legend.position = "none")

# fwrite(all_trees, "data/calibration-data/BCI/all_trees.csv")
fwrite(all_trees[census >= 1985], paste0(out_dir, "/all_trees.csv"))
# fwrite(spptable, "data/calibration-data/BCI/spptable.csv")

# merge species to all_trees
species_assigned <- unique(merge(all_trees[,.(sp, speciesID = species)], spptable[,.(sp,Genus, Species, Family)], by = "sp"))

# fwrite(species_assigned, "data/calibration-data/BCI/species_assigned.csv")
fwrite(species_assigned, paste0(out_dir,"/species_assigned.csv"))

dbh_cmTOba_m2 <- function(dbh) {
  dbh = dbh/100
  return(pi*dbh^2/4)
}



stand_dt <- all_trees[,.(
  ba = sum(dbh_cmTOba_m2(dbh_cm)*(status=="A"), na.rm = T),
  trees = sum(status=="A"),
  dbh_cm = mean(dbh_cm, na.rm = T),
  g = mean(g, na.rm = T),
  fresh_dead = sum(status=="D" & status_before=="A", na.rm = T),
  r = sum(status=="A" & status_before=="P", na.rm = T)
), by=.(census, siteID = patch, species)]

plot(table(stand_dt$r))

# define trees_before in stand_dt
stand_dt[, trees_before := shift(trees, 1, type = "lag"), by=.(siteID, species)]
stand_dt[, m := fresh_dead/trees_before, ]
stand_dt[is.infinite(m), m := 1]

# fwrite(stand_dt, "data/calibration-data/BCI/stand_dt.csv")
fwrite(stand_dt[census >= 1985], paste0(out_dir,"/stand_dt.csv"))

dbh_cm_class_size = 0.1
initial_trees <- all_trees[status=="A" & !is.na(dbh_cm), .(
  trees = sum(status == "A")
  ), by = .(dbh_cm = cut(
  dbh_cm,
  breaks = seq(1, 351, dbh_cm_class_size),
  labels = seq(1, 351-dbh_cm_class_size, dbh_cm_class_size),
  include.lowest = T), species, patch, census)
]
max(all_trees$dbh_cm, na.rm = T)
initial_trees <- initial_trees[order(patch)]
# initial_trees$cohortID = 1:nrow(initial_trees)
initial_trees[,cohortID := 1:.N,by = .(patch,census)]
initial_trees[, dbh_cm := as.numeric(as.character(dbh_cm)),]

# Npatches = uniqueN(initial_trees$patch)
#
# stand_dt[,.(trees = sum(trees, na.rm = T)/50), by = .(census)]
stand_dt[,.(ba = sum(ba, na.rm = T)/50), by = .(census)]
#
all_trees[, .(ba = sum(dbh_cmTOba_m2(dbh_cm)*(status == "A"), na.rm = T)/50), by = .(census)]
# all_trees[,.N, by = .(census, period_length)]

obs_df <- initial_trees[,.(siteID = patch, patchID = 1, cohortID, species = as.integer(as.factor(species)), dbh = dbh_cm, trees, census)]

# hist(obs_df$dbh)
obs_df[, .(ba = sum(dbh_cmTOba_m2(dbh)*trees)/50), by = .(census)]
stand_dt[, .(dbh = sum(dbh_cm*trees, na.rm = T)), by = .(census)]
obs_df[, .(dbh = sum(dbh*trees)), by = .(census)]

# stand_dt[,.(ba = sum(ba)/50), by = census]

fwrite(obs_df[census >= 1985], paste0(out_dir,"/obs_df.csv"))

