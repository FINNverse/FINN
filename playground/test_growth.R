
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Growth ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(torch)
library(FINN)
library(data.table)
library(ggplot2)
patch_size = 0.1
# trees_vec = c(1:10,10^(seq(2,4, length.out = 8)))
trees_vec = c(5)
dbh_vec = c(seq(1,100,1),seq(110,300,10))
# dbh_vec = dbh = seq(30,100, 10)

cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees = round(trees_vec*patch_size),
      dbh = dbh_vec
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)
cohort = CohortMat$new(obs_df = cohort_df1)

# dbh = cohort$dbh
# species = cohort$species
# trees = cohort$trees

sim_dt <-
  data.table(
    expand.grid(list(
      parGrowth1 = c(0.01,seq(0.2,0.8,0.2),0.99),
      parGrowth2 = seq(0,.1,.005),
      pred = seq(0,1,0.2),
      # light = seq(0,1,0.2)
      light = 0.4
    ))
  )


growth = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))) /
             (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))))

  environment = torch::torch_exp(pred) # inverse link function
  growth = shade * environment * (torch::torch_exp(-(dbh / (parGrowth[,2][species] * 100))))
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment,growth = growth) else out = growth
  return(out)
}

growth2 = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))) /
             (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))))

  environment = torch::torch_exp(pred) # inverse link function
  growth = shade * environment * (torch::torch_exp(-parGrowth[,2][species] * dbh))
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment,growth = growth) else out = growth
  return(out)
}

i=50
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  if(i==1) cat("\n")
  cat("\r", i, "of", nrow(sim_dt))
  temp_parGrowth = torch_tensor(cbind(sim_dt[i,]$parGrowth1,sim_dt[i,]$parGrowth2), requires_grad=TRUE, dtype=torch_float32())
  temp_pred = torch_tensor(array(sim_dt[i,]$pred,dim = c(1,1)), requires_grad=TRUE, dtype=torch_float32())
  light = torch_tensor(array(sim_dt[i,]$light,dim = c(1,1,1)), requires_grad=TRUE, dtype=torch_float32())
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  # m = mortality(dbh = cohort$dbh, species = cohort$species, trees = cohort$trees, parGrowth = temp_parGrowth, pred = temp_pred, light = light, debug = T)
  # dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F
  g = growth2(dbh = cohort$dbh, species = cohort$species, parGrowth = temp_parGrowth, pred = temp_pred, light = light, debug = T)
  # g = growth(dbh = cohort$dbh, species = cohort$species, parGrowth = temp_parGrowth, pred = temp_pred, light = light, debug = F)
  cohort_df1$shade = as_array(g$shade)
  cohort_df1$dbh = as_array(cohort$dbh)
  cohort_df1$environment = as.numeric(as_array(g$environment))
  cohort_df1$growth = as_array(g$growth)
  cohort_df1$pred = sim_dt[i,]$pred
  cohort_df1$parGrowth1 = sim_dt[i,]$parGrowth1
  cohort_df1$parGrowth2 = sim_dt[i,]$parGrowth2
  cohort_df1$light = sim_dt[i,]$light
  out_dt <- rbind(out_dt, cohort_df1)
  if(i==nrow(sim_dt)) cat("\n")
}

ggplot(out_dt[light == .4], aes(x = dbh, y = growth, color = factor(pred)))+
  geom_line()+
  facet_grid(paste0("sizePar=",parGrowth2)~paste0("lightPar=",parGrowth1))+
  ylim(0,NA)


ggplot(out_dt[pred == 0,
              .(shade = mean(shade)
              ), by = .(light = cut(light, seq(0, 1, 0.2), ordered_result = T, include.lowest = T), cohortID, parGrowth1)
], aes(
  x = parGrowth1,
  y = shade,
  color = light)) +
  geom_point() +
  geom_line() +
  facet_grid(light~ cohortID)

custom_palette <- c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725")

ggplot(out_dt[parMort2 == 4],
       aes(x = pred, y = light, fill = growth))+
  geom_tile()+
  facet_grid(paste0("parGrowth1=", parGrowth1)~paste0("parGrowth2=", parGrowth2))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  # scale_fill_gradientn(colors = custom_palette)
  scale_fill_gradientn(colors = custom_palette, na.value = custom_palette[length(custom_palette)])
# scale_fill_gradientn(colors = custom_palette, limits = c(0,8), na.value = custom_palette[length(custom_palette)])

ggplot(out_dt[dbh %in% c(seq(30,300,30))], aes(x = parMort2, y = gPSize, color = factor(dbh), group = dbh))+
  geom_line()+
  coord_cartesian(ylim = c(0,2))+
  geom_hline(yintercept = 0.1, col = "black", lty = 3)+
  scale_y_continuous(breaks = seq(0,2,0.1))

ggplot(out_dt[pred == 1], aes(x = dbh, y = growth, color = factor(parMort2)))+
  # geom_line()+
  geom_point()+
  facet_grid(paste0("2=",parGrowth2)~paste0("1=",parGrowth1))

ggplot(out_dt, aes(x = factor(dbh), y = growth))+
  # geom_line()+
  geom_point()+
  geom_boxplot()+
  facet_grid(paste0("2=",parGrowth2)~paste0("1=",parGrowth1))

ggplot(out_dt, aes(x = factor(light), y = growth))+
  # geom_line()+
  geom_point()+
  geom_boxplot()+
  facet_grid(paste0("2=",parGrowth2)~paste0("1=",parGrowth1))
# coord_cartesian(ylim = c(0,2))+
# geom_hline(yintercept = 0.1, col = "black", lty = 3)+
# scale_y_continuous(breaks = seq(0,2,0.1))

ggplot(out_dt, aes(x = factor(parMort1), y = shade))+
  geom_boxplot()+
  facet_wrap(~cohortID)


ggplot(out_dt, aes(
  x = parMort1,
  y = predM,
  color = pred, group = pred)) +
  geom_point() +
  # stat_smooth(se = F) +
  geom_line() +
  facet_grid(paste0("env=",environment)~paste0("light=",light))
# facet_wrap(~paste0("light=",light))

ggplot(out_dt, aes(x = factor(parMort1), y = predM))+
  geom_boxplot()

ggplot(out_dt, aes(x = factor(environment), y = predM))+
  geom_boxplot()


ggplot(out_dt[pred == 0.5 & parMort2 %in% c(1,2,3,4) & parMort1 %in% c(0,0.25,0.5,0.75, 1)], aes(x = dbh, y = light, fill = predM))+
  geom_tile()+
  facet_grid(paste0("parMort1=", parMort1)~paste0("parMort2=", parMort2))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  scale_fill_gradientn(colors = custom_palette)
# scale_fill_viridis(option = "viridis", name = "Value")
# labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "mort")


fun <- function(light, parMort1) {
  growth = pred/2 * parGrowth[,2][species] * (1-torch::torch_exp(-(dbh / (parMort[,2][species] * 100))))
}

light = seq(0,1,0.01)
plot(light,fun(light,0.5), ylim = c(0,1), type = "l", xlab = "light", ylab = "shade induced mortality")
legend("topright", title = "parMort1", legend = c("0.5","0.1","0.2","0.8","0.9"), col = c("black","red","blue","orange","green"), lty = 1)
lines(light,fun(light,.1), ylim = c(0,1), col = "red")
lines(light,fun(light,.2), ylim = c(0,1), col = "blue")
lines(light,fun(light,.8), ylim = c(0,1), col = "orange")
lines(light,fun(light,0.9), ylim = c(0,1), col = "green")
# add legend for each par and color on top left




# Define the curve function with adjustable parameters


dbh_function <- function(DBH_is, Dopt_s = 0.66, K_s = 0.14+.3, Max_s = 1) {
  DBH_is = DBH_is/10
  # Calculate the log(AGR_is + 1) based on the formula provided
  log_AGR_plus_1 <- Max_s * exp(-1/2 * (log(DBH_is / Dopt_s) / K_s)^2)

  return(log_AGR_plus_1)
}

# # Example usage of the function:
# DBH_is <- 25  # Example value for DBH_is
# Dopt_s <- 30  # Example value for Dopt_s
# K_s <- 0.5    # Example value for K_s
# Max_s <- 1.5  # Example value for Max_s
#
# result <- agr_function(DBH_is, Dopt_s, K_s, Max_s)
# print(result)

# Define the DBH values
dbh_values <- seq(0, 300, by = 1)

# Define a vector of parGrowth2 values
parGrowth2_values <- c(1, 2, 3, 4, 5)

dbh_function <- function(DBH, Dopt = 50, K = 0.7, Max = 5) {
  Max * exp(-0.5 * (log(DBH / (Dopt*10)) / K)^2)
  # torch::torch_exp(-0.5 * (log(dbh / (parGrowth[,2][species])*100) / K)^2)
}
# Plot the first curve using the first value in the parGrowth2_values vector
# plot(dbh_values, dbh_function(dbh_values, parGrowth2_values[1]), ylim = c(0, 1), type = "l", xlab = "DBH (cm)", ylab = "Density", col = "black", lty = 1)
par(mfrow = c(2,2))
for(k in c(0.1,0.5,0.8,1.5)){
  {plot(dbh_values, dbh_function(dbh_values, parGrowth2_values[1], K = k), type = "l", xlab = "DBH (cm)", ylab = "Density", col = "black", lty = 1,
        main = paste0("
       function(DBH, Dopt, K = ",k,", Max = 2)
          Max * exp(-0.5 * (log(DBH / Dopt) / K)^2)
       "))

    # Loop over the remaining parGrowth2 values and add them to the plot
    colors <- c("black", "red", "blue", "orange", "green")

    for (i in 2:length(parGrowth2_values)) {
      lines(dbh_values, dbh_function(dbh_values, parGrowth2_values[i]), col = colors[i], lty = 1)
    }
    legend("topright", title = "Dopt", legend = parGrowth2_values, col = colors, lty = 1)}
}



# Add a legend to the top right
# Define the DBH values
dbh_values <- seq(0, 400, by = 1)

# Define a vector of parGrowth2 values
parGrowth2_values <- c(30, 10, 20, 80, 200)

# Plot the first curve using the first value in the parGrowth2_values vector
plot(dbh_values, dbh_function(dbh_values, parGrowth2_values[1]), ylim = c(0, 1), type = "l", xlab = "DBH (cm)", ylab = "Density", col = "black", lty = 1)

# Loop over the remaining parGrowth2 values and add them to the plot
colors <- c("black", "red", "blue", "orange", "green")

for (i in 2:length(parGrowth2_values)) {
  lines(dbh_values, dbh_function(dbh_values, parGrowth2_values[i]), col = colors[i], lty = 1)
}

# Add a legend to the top right
legend("topright", title = "parGrowth2", legend = parGrowth2_values, col = colors, lty = 1)


growth = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))) /
             (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))))

  environment = torch::torch_exp(pred) # inverse link function
  growth = shade * environment * (torch::torch_exp(-(dbh / (parGrowth[,2][species] * 100))))
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment,growth = growth) else out = growth
  return(out)
}


