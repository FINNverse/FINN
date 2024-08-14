library(FINN)
library(data.table)
library(ggplot2)
library(torch)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## env data input ####
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



#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2. simulate data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
sp = 3L
patches = 5
sites = 10L
initCohort = CohortMat$new(dims = c(sites, patches, 10),
                           dbh = array(1, dim = c(sites, patches, 10)),
                           trees = array(1, dim = c(sites, patches, 10)),
                           sp = sp)
# torch::as_array(initCohort$dbh)
finn = FINN$new(sp = sp, env = 2L, device = "cpu", which = "all" ,
                parGrowth = matrix(c(0.1, 5), sp, 2, byrow = TRUE),
                parMort = matrix(c(runif(sp), runif(sp, 1, 4)), sp, 2, byrow = FALSE),
                parReg = runif(sp, 0.8, 0.9), # any value between 0 and 1. 0 = species needs no light for regeneration, 1 = species needs full light for regeneration
                parHeight = runif(sp, 0.3, 0.7), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                # parHeight = c(0.6,0.6,0.7), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                parGrowthEnv = list(matrix(c(10, 10, -5), sp, 2)),
                parMortEnv = list(matrix(c(-1, -1, 1)*0, sp, 2)),
                parRegEnv = list(matrix(c(1, 2, 3), sp, 2)),
                patch_size_ha = 0.1)

env = torch::torch_randn(size = c(sites, 7L, 2))*0+1

str(finn$device)

# CohortMat$new(dims = c(env$shape[1], patches, sp), sp = sp, device=finn$device)

system.time({
  FINN.seed(1)
  pred1 = finn$predict(
    # dbh = initCohort$dbh,
    # trees = initCohort$trees,
    # species = initCohort$species,
    dbh = NULL,
    trees = NULL,
    species = NULL,
    env = env, patches = patches, debug = FALSE)
})

system.time({
  FINN.seed(1)
  pred2 = finn$predict(
    # dbh = initCohort$dbh,
    # trees = initCohort$trees,
    # species = initCohort$species,
    dbh = NULL,
    trees = NULL,
    species = NULL,
    env = env, patches = patches, debug = T)
})

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 3. model output --> inventory data.frame ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# pred = pred2

results <- pred2DF(pred1, format = "long")
results <- pred2DF(pred2, format = "long")

results <- pred2DF(pred1, format = "long")
results <- pred2DF(pred2, format = "long")

str(results$site)
str(results$patch)
str(results$cohort)

ggplot(results$site[, .(value = max(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")

ggplot(results$site[, .(value = max(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 4. inventory data.frame --> model input ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

results <- pred2DF(pred2, format = "wide")



test_cohort <- CohortMat$new(results$cohort[trees > 0 & year == max(year)])

# results$cohort[trees > 0]$cohortID

test_cohort$trees

test_cohort$asDF()

dim(torch::as_array(test_cohort$trees))

dim(test_cohort$trees)
summary(test_cohort$trees)

disturbance_intensity = 0.5
disturbance_frequency = 0.02

disturbances = array(NA_integer_, dim = c(sites, timesteps, patches))
dim(disturbances)
for(site_i in 1:sites){
  for(year_i in 1:timesteps){
    dist_site_year = rbinom(n = 1, size = 1, prob = frequency)
    disturbances[site_i, year_i, ] = dist_site_year*rbinom(n = patches, size = 1, prob = intensity)
  }
}

disturbances = 1*(disturbances == 0)
disturbances_tens <- torch::torch_tensor(disturbances, device = "cpu")

test_cohort$trees*disturbances_tens[,1,]$unsqueeze(3L)


dim(disturbances_tens)

as_array(initCohort$trees)*disturbances[1, 1,] # site_i, year_i, patch_i
as_array(initCohort$trees)*disturbances[1, 1,] # site_i, year_i, patch_i

initCohort$trees*disturbances[1, 1,] # site_i, year_i, patch_i

dim(test_cohort$trees)
dim(disturbances_tens) # site_i, year_i, patch_i
dim(test_cohort$trees[,,])
dim(disturbances_tens[,1,]) # site_i, year_i, patch_i

test_cohort$trees[,,]*disturbances_tens[,1,] # site_i, year_i, patch_i


results$site

max(results$cohort$cohortID)


## 4.1 calculate response variables ####
results <- pred2DF(pred1, format = "wide")


climateDF2array(results$site, env_vars = c("dbh", "trees", "species", "ba"))



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


