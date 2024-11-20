
inputDIR <- "data/calibration-data/FIA/preparedV4/"
init_cohort1 = fread(paste0(inputDIR, "init_cohort1_dt.csv"))
species_link_dt = fread(paste0(inputDIR, "species_link_dt.csv"))

cohort1 = initCohort(init_cohort1, sp = max(species_link_dt$species))


