
library(FINN)
library(data.table)
library(ggplot2)

#' Make FINN arrays from observation data frame
#'
#' @param cohort_df a data frame with the columns Species, dbh, nTree. Optionally cohortID an be provided.
#' @param ... nothing implemented yet
#' @return list of cohort arrays dbh, nTree, Species from data frame with dendrometric data
#'
#' @export
cohort_df2arrays <- function(cohort_df){
  Nspecies = length(unique(cohort_df$Species))
  Ncohorts = nrow(cohort_df)
  Species = array(data = cohort_df$Species, dim = c(1,1,Ncohorts))
  dbh = array(data = cohort_df$dbh,dim = c(1,1,Ncohorts,1))
  nTree = array(data = cohort_df$nTree,dim = c(1,1,Ncohorts,1))
  return(list(dbh = dbh, nTree = nTree, Species = Species))
}

################################################################################
## predict with model from environment
################################################################################
species_dt = data.table(
  names = c("Fagus sylvatica", "Picea abies", "Acer pseudoplatanus", "Betula pendula", "Abies alba"),
  dbh2height = c(0.6,0.6,0.6,0.6,0.6),
  Light = c(0.1,0.15,0.5,0.6,0.07),
  Temp = c(0.6,0.3,0.7,0.4,0.3),
  Prec = c(0.6,0.7,0.4,0.4,0.5),
  maxD = c(350,300,200,100,350),
  maxG = c(3,5,7,8,4)
)

observed =fread("data/testing-data/evaluation2-obs.csv")
data = fread("data/testing-data/evaluation2-env.csv")

Nsites = 1
Ncohorts = 1

cohort_df =
  data.frame(
    cohortID = 1:5,
    Species = c(0,0,1,1,2),
    timestep = c(0,0,1,1,2),
    dbh = c(11,12,13,14,15),
    nTree = c(1,1,1,1,1)
  )

cohort_list <- cohort_df2arrays(cohort_df)

dbh2height = c(0.6,0.6,0.6)

compOut <- FINN::competition(
  dbh = cohort_list$dbh, Species = cohort_list$Species, parGlobal = dbh2height
  )


# check sensitivity of competition function

# sensitivity to dbh vs parGlobal

dbhTest = seq(0,200, length.out = 200)
parGlobalTest = seq(0.0001, 1, length.out = 9)

out_dt <- data.table()
for(dbh_i in dbhTest){
  for(parGlobal_i in parGlobalTest){
    i_cohort_df =
      data.frame(
        Species = rep(c(0),10),
        dbh = c(rep(dbh_i,5),rep(dbh_i+1,5)),
        nTree = 1
      )
    i_cohort_list = cohort_df2arrays(i_cohort_df)
    i_dbh2height = parGlobal_i
    compOut <- FINN::competition(
      dbh = i_cohort_list$dbh, Species = i_cohort_list$Species, parGlobal = i_dbh2height
    )
    out_dt <- rbind(
      out_dt,
      data.table(dbh = dbh_i,
                 parGlobal = parGlobal_i,
                 AL = as.vector(compOut),
                 heights = c(rep("small",5), rep("big", 5)))
      )

  }
}

ggplot(out_dt, aes(y = AL, x = dbh, color = heights))+
  geom_point()+
  facet_wrap(~factor(parGlobal))



str(compOut)


out_dt_all <- merge(observed, data, by = c("site", "timestep"))

circle_area <- function(d) pi * (d/100 / 2)^2
out_dt_all[,ba := (circle_area(dbh)*nTree)/0.1,]
out_dt_all[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]
out_dt_all2 <- out_dt_all[,.(
  ba = sum(ba),
  nTree = sum(nTree)
), by = .(site,timestep,SpeciesName,patch, envAvg1, envAvg2)]

p1 <- ggplot(out_dt_all2[,.(
  ba = mean(ba)
), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))+
  geom_text(aes(x = 50, y = 150, label = site), data = unique(data[,.(site, timestep = 0, envAvg1, envAvg2, SpeciesName = NA)]))
p1


m = FINN(~0 + envAct1 + envAct2,
         data = data, epochs = 1L, batch_size = 20L, start_time = 0,
         learning_rate = 1, patches = 50,
         response = observed,
         init_growth = cbind(c(0.6,0.3,0.7,0.4,0.3),
                             c(0.6,0.7,0.4,0.4,0.5)),
         init_reg = c(0.1,0.15,0.5,0.6,0.07),
         init_mort = cbind(c(0.6,0.3,0.7,0.4,0.3),
                           c(0.6,0.7,0.4,0.4,0.5))
         )
m$model$parGrowth
m$model$parReg
m$model$parMort

# model = FINN:::pkg.env$FINN$FINN(sp = 5L, env = 2L, device = "cpu")
# m2 = model$continue_fit(X = m$data$env, Y = m$data$observations, learning_rate = 0.01,
#                         batch_size = 25L, epochs = 5L, response = "ba*T")
# dim(m$data$env)
# dim(m$data$observations)

pred = FINN:::predict.FINN(m)

my_array = pred[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, data, by = c("site", "timestep"))
# out_dt_all_pred[,ba := (circle_area(dbh)*nTree)/0.1,]
out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = log(mean(ba))
), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))
p2

################################################################################
## start with correct initial cohorts
################################################################################

batch_size1 = 10
Nsites = uniqueN(data$site)
if(Nsites/batch_size1 != round(Nsites/batch_size1)) stop()

cohort = FINN:::pkg.env$FINN$CohortMat(
  dbh = array(1.0, dim = c(Nsites/batch_size1, 50, 5, 1)),
  nTree = array(1, dim = c(Nsites/batch_size1, 50, 5, 1) ),
  Species = array(matrix(0:4, 50, 5, byrow = TRUE), dim = c(Nsites/batch_size1, 50, 5))
)


m = FINN(~0 + envAct1 + envAct2,
         data = data, epochs = 1L, batch_size = 10L, start_time = 0,
         patches = 50, response = observed, cohorts = cohort)
m2 = model$continue_fit(X = m$data$env, Y = ,response = "ba*T",
                        batch_size = 100L, epochs = 1)
pred = FINN:::predict.FINN(m)

my_array = pred[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, data, by = c("site", "timestep"))
# out_dt_all_pred[,ba := (circle_area(dbh)*nTree)/0.1,]
out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = mean(ba)
), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))
p2



################################################################################
## start with correct initial cohorts
################################################################################

batch_size1 = 100
Nsites = uniqueN(data$site)
if(Nsites/batch_size1 != round(Nsites/batch_size1)) stop()

observed0 <- observed[timestep == 1]

get_array_from_obs <- function(timesteps, value, observed_dt){
  wide_dt <- dcast(observed_dt[timestep == timestep, .(patch, site,Species, dbh)], Species + patch ~ site, value.var = "dbh", fun.aggregate = sum)
  # Step 2: Convert to a matrix
  wide_matrix <- as.matrix(wide_dt[ , -c("Species", "patch"), with = FALSE])
  # Step 3: Convert to a 3D array
  # The dimensions are the unique counts of Species, patch, and site
  array_dims <- c(nrow(wide_matrix), ncol(wide_matrix), length(unique(observed_dt$Species)))
  my_array <- array(data = wide_matrix, dim = c(array_dims,1))
  return(my_array)
}

# Now observed0 is ordered by site, patch, and Species
# You can proceed with replacing the values in the array
# Assuming your 4D array is initialized as described
ntree_array <- array(1.0, dim = c(100, 50, 5, 1))

# Replace values in the array using the ordered data.table
ntree_array[cbind(observed0$site, observed0$patch, observed0$Species, 1)] <- observed0[order(site, patch, Species)]$nTree

cohort = FINN:::pkg.env$FINN$CohortMat(
  dbh = array(1.0, dim = c(100, 50, 5, 1)),
  nTree = ntree_array,
  Species = array(matrix(0:4, 50, 5, byrow = TRUE), dim = c(100, 50, 5))
)

m = FINN(~0 + envAct1 + envAct2,
         data = data, epochs = 10, batch_size = 10, start_time = 0,
         patches = 50, response = observed, CohortMat = cohort)

pred = FINN:::predict.FINN(m)

my_array = pred[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, data, by = c("site", "timestep"))
# out_dt_all_pred[,ba := (circle_area(dbh)*nTree)/0.1,]
out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = mean(ba)
), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))
p2

m2 = model$continue_fit(X = m$data$env, Y = m$data$observations, response = "ba*T",
                   batch_size = 25L, epochs = 5L)

