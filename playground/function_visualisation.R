# Load necessary library
library(ggplot2)
library(data.table)
library(torch)
library(FINN)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## create pdf ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Open the existing PDF file in append mode
# pdf("function_vis.pdf")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## height ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
test_dt =
  data.table(
    expand.grid(
      list(
        dbh = seq(0, 500,1),
        parHeight = seq(0,1,0.1)
      )
    )
  )
test_dt$rowID = 1:nrow(test_dt)

test_dt[, height := height(dbh, parHeight), by = rowID]

p_height <-
  ggplot(test_dt, aes(x = dbh, y = height, group = parHeight,color = parHeight))+
    ggtitle("height",paste0(deparse(height),collapse = "\n"))+
    ylab("height(dbh,parHeight)")+
    xlab("dbh [cm]")+
    geom_line()+
    geom_label(aes(label = parHeight, color = parHeight), test_dt[,.(height = max(height),dbh = max(dbh)),by=parHeight])
print(p_height)

#### parHeight ####
#' in principle all parameters of parHeight from 0 to 1 result in physiologicaly
#' plausible heights.
#' But very high values close to 1 are really on the edge of how tall a tree
#' can become.
#'
#' Recommendation for default parameterization.
#' The default parameterization should cover the range from 0.3 to 0.9
#'

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stem ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Generate test data
test_dt <- data.table(
  dbh = seq(0, 500, 1)
)

# Calculate basal area
test_dt[, basal_area := BA_stem(dbh)]

# Plot
p_BA_stem <- ggplot(test_dt, aes(x = dbh, y = basal_area)) +
  ggtitle("BA_stem", paste0(deparse(BA_stem), collapse = "\n")) +
  ylab("BA_stem(dbh)") +
  xlab("dbh [cm]") +
  geom_line()
print(p_BA_stem)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stand ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
trees_vec = c(0:500,10^(seq(2,4, length.out = 20)))
dbh_vec = dbh = seq(1,200, 1)
patch_size_vec = seq(0.04, 0.2, 0.02)
patch_size_vec = c(0.01,0.05,0.08,seq(0.1,1,0.2),1)

test_dt_out <- data.table()
i=patch_size_vec[1]
for(i in patch_size_vec){
  # Generate test data
  cohort_df1 =
    data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = i,
        dbh = dbh_vec
      )
    )
  cohort_df1$siteID = 1:nrow(cohort_df1)
  cohort_df1$trees = round(cohort_df1$trees_ha*i)
  cohort = CohortMat$new(obs_df = cohort_df1)

  basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = i)
  light = competition(dbh = cohort$dbh, species = cohort$species,
                      trees = cohort$trees, parHeight = torch::torch_tensor(c(0.5)),
                      patch_size_ha = i,
                      h = 0)
  cohort_df1$basal_area = torch::as_array(basal_area)[,1,1]
  cohort_df1$light = torch::as_array(light)[,1,1]
  cohort_df1$patch_size_ha = i
  test_dt_out <- rbind(test_dt_out,cohort_df1)
}

max_ba = 60
test_dt_out[basal_area > max_ba, basal_area := max_ba,]
p_BA_stand <- ggplot(test_dt_out[trees_ha <=500],aes(y = factor(trees_ha), x = dbh, fill = basal_area))+
  geom_tile()+
  ylab("trees [N/ha]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="basal area\n[m^2/ha]"))+
  ggtitle("BA_stand",paste0(deparse(BA_stand),collapse = "\n"))+
  scale_fill_viridis_c(breaks = seq(0,max_ba,10), labels = c(seq(0,max_ba-10,10),paste0(">",max_ba)))+
  scale_y_discrete(breaks = seq(0,500,100))+
  facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ",x, " ha")))
print(p_BA_stand)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## competition ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# sensitivity of competition to parHeight (species)
parHeight_vec = c(0.3,0.7)
for(i in parHeight_vec){
  test_dt_out <- data.frame()
  i=patch_size_vec[1]
  # Generate test data
  cohort_df1 =
    data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = 0,
        dbh = dbh_vec
      )
    )
  cohort_df1$siteID = 1:nrow(cohort_df1)
  cohort_df1$trees = round(cohort_df1$trees_ha*i)
  cohort = CohortMat$new(obs_df = cohort_df1)

  }

# sensitivity of competition to h: light availability vs. height

trees_vec = c(0:500,10^(seq(2,4, length.out = 20)))
dbh_vec = dbh = seq(1,200, 1)
patch_size_vec = seq(0.04, 0.2, 0.02)
patch_size_vec = c(0.1)



test_dt_out <- data.table()
# Generate test data
cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees_ha = trees_vec,
      patch_size_ha = 0.1,
      dbh = dbh_vec
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)
cohort_df1$trees = round(cohort_df1$trees_ha*i)
cohort = CohortMat$new(obs_df = cohort_df1)

basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = 0.1)
light = competition(dbh = cohort$dbh, species = cohort$species,
                    trees = cohort$trees, parHeight = torch::torch_tensor(parHeight_vec),
                    h = 1, patch_size_ha = 0.1)
cohort_df1$basal_area = torch::as_array(basal_area)[,1,1]
cohort_df1$light = torch::as_array(light)[,1,1]
test_dt_out <- rbind(test_dt_out,cohort_df1)

p_competition1 <- ggplot(test_dt_out[trees_ha <=500],aes(y = factor(trees_ha), x = dbh, fill = light))+
  geom_tile()+
  ylab("trees [N/ha]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light[%]"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))+
  # scale_fill_viridis_c(breaks = seq(0,max_ba,10), labels = c(seq(0,max_ba-10,10),paste0(">",max_ba)))+
  scale_y_discrete(breaks = seq(0,500,100))+
  facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ",x, " ha")))
print(p_competition1)


patch_size = 0.1

trees_vec = c(1:99,10^(seq(2,4, length.out = 16)))
dbh_vec = dbh = seq(1,300, 1)

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

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

light = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0, patch_size_ha = 0.1)
cohort_df1$light = torch::as_array(light)[,1,1]

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light*100))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## regeneration ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

patch_size = 0.1
# trees_vec = c(1:10,10^(seq(2,4, length.out = 8)))
trees_vec = c(0:100)
dbh_vec = 80
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

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

light = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0, patch_size_ha = 0.1)
cohort_df1$light = torch::as_array(light)[,1,1]
hist(cohort_df1$light)

sim_dt <-
  data.table(
    expand.grid(list(
      parReg = seq(0,1,0.02),
      pred = seq(0,500,5)
    ))
  )

# r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
# regeneration2 = function(species, parReg, pred, light) {
#   if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
#   environment = pred
#   regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
#   mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2 )
#   regeneration1 = sample_poisson_relaxed(mean, num_samples = 1000) # TODO, check if exp or not?! lambda should be always positive!
#   regeneration2 = regeneration1 + regeneration1$round()$detach() - regeneration1$detach()
#   return(list(regP = regP, mean = mean, regeneration1 = regeneration1, regeneration2 = regeneration2))
# }

light1 = 0.5
parReg = seq(0,1,0.001)
plot(parReg,plogis((light1 + (1-parReg) - 1)/1e-1))
plot(parReg,(light1 + (1-parReg) - 1), ylim = c(0,1))

ffun <- function(light, parMort1) {
  # light=1-light
  # parMort1 = 1-parMort1
  # 1-(1-(1 / (1 + exp(-10 * (light - parMort1))) - 1 / (1 + exp(10 * parMort1))) / (1 - 1 / (1 + exp(10 * (1 - parMort1)))))
  # 1 - ((1 / (1 + exp(-10 * (light - parMort1))) - 1 / (1 + exp(10 * parMort1))) / (1 / (1 + exp(-10 * (1 - parMort1))) - 1 / (1 + exp(10 * parMort1))))
  1 - ((1 / (1 + exp(-20 * (light - parMort1))) - 1 / (1 + exp(20 * parMort1))) / (1 / (1 + exp(-20 * (1 - parMort1))) - 1 / (1 + exp(20 * parMort1))))
}
ffun <- function(light, parMort1, steepness = 100) {
  1 - ((1 / (1 + exp(-steepness * (light - parMort1))) - 1 / (1 + exp(steepness * parMort1))) /
         (1 / (1 + exp(-steepness * (1 - parMort1))) - 1 / (1 + exp(steepness * parMort1))))
}

ffun <- function(light, parMort1, base_steepness = 10) {
  # Scale steepness towards the edges
  scaled_steepness <- base_steepness / (0.5 - abs(parMort1 - 0.5))

  1 - ((1 / (1 + exp(-scaled_steepness * (light - parMort1))) - 1 / (1 + exp(scaled_steepness * parMort1))) /
         (1 / (1 + exp(-scaled_steepness * (1 - parMort1))) - 1 / (1 + exp(scaled_steepness * parMort1))))
}

steep_plot = 10

# Test the function
# Sorted parReg values
sorted_parReg <- c(0.1, 0.2, 0.4, 0.5, 0.9)
sorted_colors <- c("red", "blue", "green", "black", "orange")  # Corresponding sorted colors

# Test the function with sorted values
light <- seq(0, 1, 0.01)
plot(light, ffun(light, sorted_parReg[4], steep_plot), ylim = c(0, 1), type = "l", xlab = "light", ylab = "output", col = sorted_colors[4])
lines(light, ffun(light, sorted_parReg[1], steep_plot), col = sorted_colors[1])
lines(light, ffun(light, sorted_parReg[2], steep_plot), col = sorted_colors[2])
lines(light, ffun(light, sorted_parReg[5], steep_plot), col = sorted_colors[5])
lines(light, ffun(light, sorted_parReg[3], steep_plot), col = sorted_colors[3])

# Add legend for each par and color on top left
legend("bottomleft", title = "parMort1", legend = sorted_parReg, col = sorted_colors, lty = 1)


light = seq(0,1,0.01)
plot(light,fun(light,0.5), ylim = c(0,1), type = "l", xlab = "light")
lines(light,fun(light,.1), col = "red")
lines(light,fun(light,.2), col = "blue")
lines(light,fun(light,.9), col = "orange")
lines(light,fun(light,0.4), col = "green")
# add legend for each par and color on top left
legend("topleft", title = "parReg", legend = c("0.5","0.1","0.2","0.9","0.4"), col = c("black","red","blue","orange","green"), lty = 1)

# legend("topleft", title = "parReg", legend = c("0.5","0.1","0.2","0.8","0.9"), col = c("black","red","blue","orange","green"), lty = 1)



plot(x,smooth_transition(x,0.9), col = "green")
smooth_transition(x,.001)

plot(x,exp(1/x))

i=1
library(torch)
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = light, patch_size_ha = 0.1, debug = T)
  cohort_df1$regP = as_array(r$regP)
  cohort_df1$mean = as_array(r$mean)
  cohort_df1$r1 = as_array(r$regeneration1)
  cohort_df1$r2 = as_array(r$regeneration2)
  cohort_df1$parReg = sim_dt[i,]$parReg
  cohort_df1$pred = sim_dt[i,]$pred
  out_dt <- rbind(out_dt, cohort_df1)
}

ggplot(out_dt, aes(x = pred, y = parReg, fill = mean))+
  geom_tile()+
  facet_wrap(~paste0("light=", cut(light,6)))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "lambda (mean of poisson)")

ggplot(out_dt, aes(x = pred, y = parReg, fill = r2))+
  geom_tile()+
  facet_wrap(~paste0("light=", cut(light,6)))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "simulated regeneration per patch")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## mortality ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(torch)
library(FINN)
library(data.table)
library(ggplot2)
patch_size = 0.1
# trees_vec = c(1:10,10^(seq(2,4, length.out = 8)))
trees_vec = c(50)
dbh_vec = seq(10,300,10)
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
#
# cohort_df1b =
#   data.frame(
#     patchID = 1,
#     cohortID = 2,
#     species = 1,
#     expand.grid(
#       trees = round(trees_vec*patch_size),
#       dbh = dbh_vec+1
#     )
#   )
# cohort_df1b$siteID = 1:nrow(cohort_df1b)
#
# cohort_df1 <- data.table(rbind(cohort_df1, cohort_df1b))

cohort = CohortMat$new(obs_df = cohort_df1)

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

# light1 = competition(cohort$dbh, cohort$species, cohort$trees,
#                     parHeight = torch::torch_tensor(0.5), patch_size_ha = 0.1)
# dim(torch::as_array(light1))
# cohort_df1[["light1"]] = c(torch::as_array(light1)[,1,1],torch::as_array(light1)[,1,2])
# cohort_df1$light2 = torch::as_array(light1)[,1,2]
# hist(cohort_df1$light1)
# hist(cohort_df1$light2)
# hist(cohort_df1[cohortID == 1]$light1)
# hist(cohort_df1$light2)
# plot(cohort_df1$light1, cohort_df1$light2, xlim = range(as_array(light1)), ylim = range(as_array(light1)))
#
# light2 = competition(cohort$dbh, cohort$species, cohort$trees,
#                     parHeight = torch::torch_tensor(0.5), patch_size_ha = 0.1)
# dim(torch::as_array(light1))
# dim(torch::as_array(light2))
# cohort_df1$light2 = torch::as_array(light2)[,1,2]
# hist(cohort_df1$light2)
# hist(cohort_df1$light)

# gPSize = torch_clamp(0.1*(dbh/torch_clamp((parMort[,2][species]*100), min = 0.00001))$pow(2.3), max = 1.0)
# 1-torch::torch_exp(-(dbh / (parMort[,2][species] * 100))^1)
# # calculate gPsizes without torch_clamp
# gPSize = (dbh/(parMort[,2][species]*100)^(2.3)
#
#
#
dbh = seq(0,400,1)
parMort2 = seq(0,.4,0.1)
fun <- function(dbh,parMort2) (dbh/(parMort2*0.1*100))^(2.3)

fun <- function(dbh, parMort2) {
  (exp(-(dbh / (parMort2 * 100))^1))
}

plot(dbh,fun(dbh,1), ylim = c(0,2), type = "l", xlab = "dbh", ylab = "gPSize")
lines(dbh,fun(dbh,0.1), col = "red")
lines(dbh,fun(dbh,0.5), col = "blue")
lines(dbh,fun(dbh,1), col = "green")
lines(dbh,fun(dbh,2), col = "orange")
lines(dbh,fun(dbh,3), col = "black")
lines(dbh,fun(dbh,4), col = "purple")
abline(h = 1, col = "black", lty = 3)
# add legend for each par and color on top left
legend("topright", title = "parMort2", legend = c("1","0.1","0.5","2","3","4"), col = c("black","red","blue","orange","green","purple"), lty = 1)


sim_dt <-
  data.table(
    expand.grid(list(
      parMort1 = c(0.01,seq(0.2,0.8,0.2),0.99),
      parMort2 = seq(0,4,1),
      pred = seq(0,1,0.1),
      light = seq(0,1,0.1)
    ))
  )


i=50
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  if(i==1) cat("\n")
  cat("\r", i, "of", nrow(sim_dt))
  temp_parMort = torch_tensor(cbind(sim_dt[i,]$parMort1,sim_dt[i,]$parMort2), requires_grad=TRUE, dtype=torch_float32())
  temp_pred = torch_tensor(array(sim_dt[i,]$pred,dim = c(length(as_array(cohort$dbh)),1)), requires_grad=TRUE, dtype=torch_float32())
  light = torch_tensor(array(sim_dt[i,]$light,dim = c(length(as_array(cohort$dbh)),1,1)), requires_grad=TRUE, dtype=torch_float32())
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  m = mortality(dbh = cohort$dbh, species = cohort$species, trees = cohort$trees, parMort = temp_parMort, pred = temp_pred, light = light, debug = T)
  cohort_df1$shade = as_array(m$shade)
  cohort_df1$dbh = as_array(cohort$dbh)
  cohort_df1$environment = as_array(m$environment)
  cohort_df1$gPSize = as_array(m$gPSize)
  cohort_df1$predM = as_array(m$predM)
  cohort_df1$mort1 = as_array(m$mort1)
  cohort_df1$mort2 = as_array(m$mort2)
  cohort_df1$pred = sim_dt[i,]$pred
  cohort_df1$parMort1 = sim_dt[i,]$parMort1
  cohort_df1$parMort2 = sim_dt[i,]$parMort2
  cohort_df1$light = sim_dt[i,]$light
  out_dt <- rbind(out_dt, cohort_df1)
  if(i==nrow(sim_dt)) cat("\n")
}

ggplot(out_dt[pred == 0,
              .(shade = mean(shade)
              ), by = .(light = cut(light, seq(0, 1, 0.2), ordered_result = T, include.lowest = T), cohortID, parMort1)
              ], aes(
  x = parMort1,
  y = shade,
  color = light)) +
  geom_point() +
  geom_line() +
  facet_grid(light~ cohortID)

custom_palette <- c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725")

ggplot(out_dt[pred == 10 & parMort2 %in% c(1,2,3,4) & parMort1 %in% c(0,0.25,0.5,0.75, 1)],
       aes(x = dbh, y = light, fill = predM))+
  geom_tile()+
  facet_grid(paste0("parMort1=", parMort1)~paste0("parMort2=", parMort2))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  # scale_fill_gradientn(colors = custom_palette)
  scale_fill_gradientn(colors = custom_palette, limits = c(0,1), na.value = custom_palette[length(custom_palette)])

summary(out_dt$predM)
ggplot(out_dt,
       aes(x = pred, y = light, fill = predM))+
  geom_tile()+
  facet_grid(paste0("parMort1=", parMort1)~paste0("parMort2=", parMort2))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  # scale_fill_gradientn(colors = custom_palette)
  scale_fill_gradientn(colors = custom_palette, limits = c(0,.1), na.value = custom_palette[length(custom_palette)])

ggplot(out_dt[dbh %in% c(seq(30,300,30))], aes(x = parMort2, y = gPSize, color = factor(dbh), group = dbh))+
  geom_line()+
  coord_cartesian(ylim = c(0,2))+
  geom_hline(yintercept = 0.1, col = "black", lty = 3)+
  scale_y_continuous(breaks = seq(0,2,0.1))

ggplot(out_dt[dbh %in% c(seq(30,300,30))], aes(x = dbh, y = gPSize, color = factor(parMort2)))+
  geom_line()+
  coord_cartesian(ylim = c(0,2))+
  geom_hline(yintercept = 0.1, col = "black", lty = 3)+
  scale_y_continuous(breaks = seq(0,2,0.1))

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
  light = 1-light
  ((((1 / (1 + torch_exp(-10 * (light - parMort1))) - 1 / (1 + torch_exp(10 * parMort1))) / (1 - 1 / (1 + torch_exp(10 * (1 - parMort1)))))) * (1 - light) + light)
}

light = seq(0,1,0.01)
plot(light,fun(light,0.5), ylim = c(0,1), type = "l", xlab = "light", ylab = "shade induced mortality")
legend("topright", title = "parMort1", legend = c("0.5","0.1","0.2","0.8","0.9"), col = c("black","red","blue","orange","green"), lty = 1)
lines(light,fun(light,.1), ylim = c(0,1), col = "red")
lines(light,fun(light,.2), ylim = c(0,1), col = "blue")
lines(light,fun(light,.8), ylim = c(0,1), col = "orange")
lines(light,fun(light,0.9), ylim = c(0,1), col = "green")
# add legend for each par and color on top left

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
dbh_vec = seq(10,300,10)
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

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

sim_dt <-
  data.table(
    expand.grid(list(
      parGrowth1 = c(0.01,seq(0.2,0.8,0.2),0.99),
      parGrowth2 = seq(1,8,1),
      parMort1 = 0,
      parMort2 = seq(0,4,1),
      pred = seq(0,1,0.2),
      light = seq(0,1,0.2)
    ))
  )


i=50
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  if(i==1) cat("\n")
  cat("\r", i, "of", nrow(sim_dt))
  temp_parMort = torch_tensor(cbind(sim_dt[i,]$parMort1,sim_dt[i,]$parMort2), requires_grad=TRUE, dtype=torch_float32())
  temp_parGrowth = torch_tensor(cbind(sim_dt[i,]$parGrowth1,sim_dt[i,]$parGrowth2), requires_grad=TRUE, dtype=torch_float32())
  temp_pred = torch_tensor(array(sim_dt[i,]$pred,dim = c(length(as_array(cohort$dbh)),1)), requires_grad=TRUE, dtype=torch_float32())
  light = torch_tensor(array(sim_dt[i,]$light,dim = c(length(as_array(cohort$dbh)),1,1)), requires_grad=TRUE, dtype=torch_float32())
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  # m = mortality(dbh = cohort$dbh, species = cohort$species, trees = cohort$trees, parGrowth = temp_parGrowth, pred = temp_pred, light = light, debug = T)

  g = growth(dbh = cohort$dbh, species = cohort$species, parGrowth = temp_parGrowth, parMort = temp_parMort, pred = temp_pred, light = light, debug = T)
  cohort_df1$shade = as_array(g$shade)
  cohort_df1$dbh = as_array(cohort$dbh)
  cohort_df1$environment = as_array(g$environment)
  cohort_df1$growth = as_array(g$growth)
  cohort_df1$pred = sim_dt[i,]$pred
  cohort_df1$parGrowth1 = sim_dt[i,]$parGrowth1
  cohort_df1$parGrowth2 = sim_dt[i,]$parGrowth2
  cohort_df1$parMort2 = sim_dt[i,]$parMort2
  cohort_df1$light = sim_dt[i,]$light
  out_dt <- rbind(out_dt, cohort_df1)
  if(i==nrow(sim_dt)) cat("\n")
}

ggplot(out_dt[pred == 0,
              .(shade = mean(shade)
              ), by = .(light = cut(light, seq(0, 1, 0.2), ordered_result = T, include.lowest = T), cohortID, parMort1)
], aes(
  x = parMort1,
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
