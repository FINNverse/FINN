library(FINN)
library(data.table)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 1. env data input ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Define the number of days in each month
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Define the years and site IDs
years <- 1970:2019
site_ids <- 1:4

# Create an empty list to store data
climate_data_list <- list()

# Loop over each year, month, and site to generate data
for (site_id in site_ids) {
  for (year in years) {
    for (month in 1:12) {
      days <- 1:days_in_month[month]

      # Handle leap years
      if (month == 2 && (year %% 4 == 0 && (year %% 100 != 0 || year %% 400 == 0))) {
        days <- 1:29
      }

      # Create a data frame for the current site, year, and month
      climate_data <- data.frame(
        siteID = site_id,
        year = year,
        month = month,
        day = days,
        tmp = runif(length(days), -10, 30),  # Random temperature values
        pre = runif(length(days), 0, 200)    # Random precipitation values
      )

      # Add to the list
      climate_data_list <- append(climate_data_list, list(climate_data))
    }
  }
}

# Combine the list into a single data frame
climate_dt <- rbindlist(climate_data_list)

# Use the function to create the array
climate_dt_day <- climate_dt
resultDay <- climateDF2array(climate_dt_day, env_vars = c("tmp", "pre"))
str(resultDay)

climate_dt_month <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID,year,month)]
resultMonth <- climateDF2array(climate_dt = data.frame(climate_dt_month), env_vars = c("tmp", "pre"))
str(resultMonth)

climate_dt_year <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID,year)]
resultYear <- climateDF2array(climate_dt = climate_dt_year, env_vars = c("tmp", "pre"))
str(resultYear)

# Create arrays
resultDay <- climateDF2array(climate_dt, env_vars = c("tmp", "pre"))
resultMonth <- climateDF2array(climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID, year, month)], env_vars = c("tmp", "pre"))
resultYear <- climateDF2array(climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID, year)], env_vars = c("tmp", "pre"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2. simulate data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
sp = 3L
patches = 40L
sites = 20L
initCohort = CohortMat$new(dims = c(sites, patches, 10),
                           dbh = array(1, dim = c(sites, patches, 10)),
                           trees = array(1, dim = c(sites, patches, 10)),
                           sp = sp)


finn = FINN$new(sp = sp, env = 2L, device = "cpu", which = "all" ,
                parGrowth = matrix(c(0.1, 5), sp, 2, byrow = TRUE),
                parMort = matrix(c(runif(sp), runif(sp, 1, 4)), sp, 2, byrow = FALSE),
                parReg = runif(sp, 0.8, 0.9), # any value between 0 and 1. 0 = species needs no light for regeneration, 1 = species needs full light for regeneration
                # parHeight = runif(sp, 0.3, 0.7), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                parHeight = c(0.6,0.6,0.7), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                parGrowthEnv = list(matrix(c(10, 10, -5), sp, 2)),
                parMortEnv = list(matrix(c(-1, -1, 1)*0, sp, 2)),
                parRegEnv = list(matrix(c(1, 2, 3), sp, 2)),
                patch_size_ha = 0.1)
env = torch::torch_randn(size = c(sites, 2000L, 2))
# env = torch::torch_zeros(size = c(sites, 100, 2))

system.time({
  pred = finn$predict(dbh = initCohort$dbh,
                      trees = initCohort$trees,
                      species = initCohort$species,
                      response = "BA*T", env = env, patches = patches, debug = FALSE)
})

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 3. model output --> inventory data.frame ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# 1 -> dbh/ba, 2 -> counts, 3 -> AL, 4 -> growth rates, 5 -> mort rates, 6 -> reg rates

# Example structure of `pred` and `out_names`
# Assuming pred is already defined as given in the question.
out_names <- c("dbh_ba", "trees", "AL", "growth", "mort", "reg")

# Initialize an empty data.table to store the final result
inventory_dt <- NULL

# Loop through pred and convert each to a data frame using array2DF
for (i in seq_along(pred)) {
  i_name <- out_names[i]

  # Convert array to a data frame using array2DF
  df <- data.table(
    array2DF(
      torch::as_array(pred[[i]]),
      responseName = i_name, base = c("siteID", "year", "species"),
      simplify = F, allowLong = T)
    )

  # Rename the dimension columns to match siteID, year, species
  setnames(df, old = c("Var1", "Var2", "Var3"), new = c("siteID", "year", "species"))

  # Convert siteID, year, species to numeric (or integer)
  df[, siteID := as.integer(gsub("siteID","",siteID))+1]
  df[is.na(siteID), siteID := 1]
  df[, year := as.integer(gsub("year","",year))+1]
  df[is.na(year), year := 1]
  df[, species := as.integer(gsub("species","",species))+1]
  df[is.na(species), species := 1]

  # If inventory_dt is NULL, initialize it with the first data frame
  if (is.null(inventory_dt)) {
    inventory_dt <- df
  } else {
    # Otherwise, merge the new data frame with the existing inventory_dt
    inventory_dt <- merge(inventory_dt, df, by = c("siteID", "year", "species"))
  }
}

# Show the structure of inventory_dt to confirm
summary(inventory_dt)

melt_dt <- melt(inventory_dt, id.vars = c("siteID","year", "species"), measure.vars = out_names)

ggplot(melt_dt[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = out_names[i]) +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 4. inventory data.frame --> model input ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

## 4.1 calculate response variables ####

# growth rates

# mortality rates

# regeneration rates

# basal area

# N trees

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 5. fit model  ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 6. true vs. fitted  ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=


