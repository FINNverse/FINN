library(data.table)

census_years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)

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

fwrite(env_dt, "data/calibration-data/BCI/env_dt.csv")

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

all_trees <- all_trees[order(census)]
all_trees[grepl("A",status), status := "A",]
all_trees[,species := as.integer(as.factor(sp)),]
all_trees[, ":="(
  status_before = data.table::shift(status, 1, type = "lag"),
  dbh_before = data.table::shift(dbh, 1, type = "lag")
  ),
  by = .(stemID)]

all_trees <- all_trees[!is.na(gx) & !is.na(gy)]
all_trees[,":="(
  x_class = cut(gx, breaks = seq(0, 1000, 100), labels = LETTERS[1:length(seq(0, 900, 100))], include.lowest = T),
  y_class = cut(gy, breaks = seq(0, 500, 100), labels = LETTERS[1:length(seq(0, 400, 100))], include.lowest = T)
),]

all_trees[, patch := as.integer(as.factor(paste0(x_class,y_class))),]

stand_dt <- all_trees[,.(
  ba = round(sum(ba*(status=="A")),4),
  trees = sum(status=="A"),
  dbh = mean(dbh/10),
  g = mean(dbh/10-dbh_before/10, na.rm = T),
  fresh_dead = sum(status=="D" & status_before=="A", na.rm = T),
  r = sum(status=="A" & status_before=="P", na.rm = T)
), by=.(census, siteID = patch, species)]

# define trees_before in stand_dt
stand_dt[, trees_before := shift(trees, 1, type = "lag"), by=.(siteID, species)]
stand_dt[, m := fresh_dead/trees_before, ]
stand_dt[is.infinite(m), m := 1]

fwrite(stand_dt, "data/calibration-data/BCI/stand_dt.csv")

dbh_class_size = 1
initial_trees <- all_trees[status=="A" & !is.na(dbh) & census == 1982, .(
  trees = sum(status == "A")
  ), by = .(dbh = cut(
  dbh/10,
  breaks = seq(1, 350, dbh_class_size),
  labels = seq(1, 349, dbh_class_size),
  include.lowest = T), species, patch)
]

initial_trees <- initial_trees[order(patch)]
initial_trees$cohortID = 1:nrow(initial_trees)

initial_trees[, dbh := as.numeric(as.character(dbh)),]

obs_df <- initial_trees[,.(siteID = patch, patchID = 1, cohortID, species = as.integer(as.factor(species)), dbh, trees)]

fwrite(obs_df, "data/calibration-data/BCI/obs_df.csv")
obs_array <- FINN::obsDF2arrays(obs_df)
cohort1 <- FINN::CohortMat$new(obs_df)


