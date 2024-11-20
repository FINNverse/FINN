library(data.table)
library(FINN)
library(torch)
library(forcats)
library(dplyr)
library(ggplot2)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## read and format parameters ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

pars <- readRDS("data/calibration-output/BCI_parameters_01ha_9_9_wo_T_v2.RDS")
env_dt = fread("data/calibration-data/BCI-1h-patch-V2/env_dt.csv")
env_dt = env_dt[,-"T_mean"]

pars_i = pars[[length(pars)]]
env_names = c("Inter", colnames(env_dt)[-c(1)])
out <- data.table()
for(i in 1:length(pars_i)){
  par_name = names(pars_i)[i]
  # print(i)
  # print(par_name)
  dt <- data.table(pars_i[[i]])
  col_names = paste0(par_name, "__", 1:ncol(pars_i[[i]]))
  if(grepl("nn", par_name)){
    col_names = paste0(par_name, "__", env_names)
  }
  if(grepl("parHeight", par_name)) {
    col_names = paste0(par_name)
  }
  if(grepl("parMort", par_name) | grepl("parGrowth", par_name)) {
    col_names = paste0(par_name, "__", c("light", "size"))
  }
  if(grepl("parReg", par_name)) {
    col_names = paste0(par_name, "__", c("light"))
  }
  colnames(dt) <- col_names
  out = cbind(out, dt)
}

out$speciesID = 1:nrow(out)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## PCA analysis ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=


cor(out[,grep("par|Inter", colnames(out)), with = FALSE])

# PCA all species
data = out[,grep("par|Inter", colnames(out)), with = FALSE]
names(data)[names(data) == "nnReg.0.weight__Inter"] <- "Regeneration (Intercept)"
names(data)[names(data) == "nnMort.0.weight__Inter"] <- "Mortality (Intercept)"
names(data)[names(data) == "nnGrowth.0.weight__Inter"] <- "Growth (Intercept)"
# data = out[,grep("par", colnames(out)), with = FALSE]

# Perform PCA on the data
pca_result <- prcomp(data, scale. = TRUE)
summary(pca_result)

#show explained variance in % for PC1 and PC2
explained_variance = pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

set.seed(123)
# kmeans_result <- kmeans(pca_result$x[, 1:2], centers = 5, nstart = 25)
kmeans_result <- kmeans(pca_result$x[, 1:2],
                        centers = rbind(
                          c(1,0),
                          c(-1,0),
                          c(0,-1),
                          c(0,1),
                          c(0,0)
                          ), nstart = 25)

# Create a data frame with the PCA results and cluster assignments
pca_data <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2],
                       cluster = as.factor(kmeans_result$cluster)
)

# Extract loadings for the first two principal components
loadings <- pca_result$rotation[, 1:2]
loadings_df <- data.frame(Variable = rownames(loadings),
                          PC1 = loadings[, 1] * 5,  # Scaling for better visualization
                          PC2 = loadings[, 2] * 5)


colors <- c(
  "Slow" = "#E07A9E",        # Pink
  "Fast" = "#F19D38",        # Orange
  "Long-lived pioneer (LLP)" = "#5E9F6C",  # Green
  "Short-lived breeder (SLB)" = "#4F7CBF", # Blue
  "Intermediate" = "#DF7045"  # Orange-red (close to red)
)

pca_data$PFT = factor(
  pca_data$cluster,
  levels = c(1,2,3,4,5),
  labels = c("Slow","Fast","Short-lived breeder (SLB)","Long-lived pioneer (LLP)","Intermediate")
  )

loadings_df$Variable2 =
  factor(
    loadings_df$Variable,
    levels = c("parHeight","parGrowth__light", "parGrowth__size", "parMort__light", "parMort__size", "parReg__light", "Regeneration (Intercept)"  ,  "Growth (Intercept)", "Mortality (Intercept)"),
    labels = c("Height", "Growth (light)", "Growth (size)", "Mortality (light)", "Mortality (size)", "Regeneration (light)", "Regeneration (Intercept)"   , "Growth (Intercept)", "Mortality (Intercept)")
  )
# levels(loadings_df$Variable2)[levels(loadings_df$Variable2) == "nnReg.0.weight__Inter"] <- "Regeneration (Intercept)"
# levels(loadings_df$Variable2)[levels(loadings_df$Variable2) == "nnMort.0.weight__Inter"] <- "Mortality (Intercept)"
# levels(loadings_df$Variable2)[levels(loadings_df$Variable2) == "nnGrowth.0.weight__Inter"] <- "Growth (Intercept)"
# Create the plot
p = ggplot(pca_data, aes(x = PC1, y = PC2, color = PFT)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey")+
  geom_point(size = 3, alpha = 1) +
  theme_minimal() +
  labs(
    x = paste0("Growth-survival tradeoff (", round(explained_variance[1]),"%)"),
    y = paste0("Stature-recruitment tradeoff (", round(explained_variance[2]),"%)")
    ) +
  theme(legend.position = "right") +
  # scale_color_manual(values = c("magenta", "cyan", "green", "blue", "orange"),
  #                    labels = c("Cluster 1", "Cluster 2", "Cluster 3",
  #                               "Cluster 4", "Cluster 5")) +
  scale_color_manual(values = colors,
                     labels = names(colors)) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.4, "cm")), color = "black", linewidth = 1.5) +
  ggrepel::geom_text_repel(
    data = loadings_df[loadings_df$PC1 < 0, ],
    aes(x = PC1, y = PC2, label = Variable2), min.segment.length = 10,
    hjust = 1, vjust = -1, color = "black", fontface = "bold")+
  ggrepel::geom_text_repel(
    data = loadings_df[loadings_df$PC1 >= 0, ],
    aes(x = PC1, y = PC2, label = Variable2), , min.segment.length = 10,
    hjust = 0, vjust = 0.3, color = "black", fontface = "bold")+
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-3, 4))+
  ggthemes::theme_base()
# geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = Variable),
#           hjust = 0.5, vjust = 0.5, color = "black")
p
# save p as svg
ggsave("gfoe2024/BCI_PCA_plot.png", p, width = 23, height = 15, units = "cm")

assigned_species <- fread("data/calibration-data/BCI/species_assigned.csv")
ruger_species <- data.table(readxl::read_xlsx("data/calibration-output/aaz4797_ruger_data_s1.xlsx", sheet = 2))

pca_data$speciesID = as.integer(rownames(pca_data))

pca_data <- merge(pca_data, assigned_species, by = "speciesID")
ruger_species$sp <- tolower(ruger_species$sp)

pca_data2 <- data.table(merge(pca_data, ruger_species[,.(sp, PFT_1axis, PFT_2axes, PC1score, PC2score)], by = "sp", all = T))

par(mfrow = c(1, 2))
model_dat <- pca_data2[!is.na(PC1score)][order(PC1score)]
plot(PC1 ~ PC1score, model_dat, xlab = "Rüger et al. (2020) PC1", ylab = "FINN PC1")
fit <- lm(PC1 ~ PC1score, model_dat)
pred <- predict(fit, newdata = model_dat, interval = "confidence")
lines(model_dat$PC1score, pred[,1], col = "red")
# # add confidence interval
lines(model_dat$PC1score, pred[,2], col = "red", lty = 2)
lines(model_dat$PC1score, pred[,3], col = "red", lty = 2)
summary(fit)

# same for PC2
model_dat2 <- pca_data2[!is.na(PC2score)][order(PC2score)]
plot(PC2 ~ PC2score, model_dat2, xlab = "Rüger PC2", ylab = "FINN PC2")
fit2 <- lm(PC2 ~ PC2score, model_dat2)
pred2 <- predict(fit2, newdata = model_dat2, interval = "confidence")
lines(model_dat2$PC2score, pred2[,1], col = "red")
# add confidence interval
lines(model_dat2$PC2score, pred2[,2], col = "red", lty = 2)
lines(model_dat2$PC2score, pred2[,3], col = "red", lty = 2)
summary(fit2)







create_plot = function(data, env = "", group = tree_families) {

  pp = pars[[length(pars)]]
  R_growth = data[,-c(1)]**2 / rowSums(data[,-1]**2)
  df = data.frame(R2 = c(R_growth[,1], R_growth[,2], R_growth[,3], R_growth[,4], R_growth[,5]),
                  part = rep(c("Prec", "SR_kW_m2", "RH_prc", "T_max", "T_min"), each = nrow(R_growth)),
                  group = rep(group, 5)
  )

  df_normalized <- df %>%
    group_by(group) %>%
    mutate(R2_norm = R2 / sum(R2))

  df_normalized$env = env
  return(df_normalized)
}

pp = pars[[length(pars)]]
results =
  rbind(create_plot(pp$nnMort.0.weight, env = "Mortality", group = kmeans_result$cluster),
        create_plot(pp$nnGrowth.0.weight, env = "Growth", group = kmeans_result$cluster),
        create_plot(pp$nnReg.0.weight, env = "Regeneration", group = kmeans_result$cluster))

# Normalize R2 values to sum up to 1 within each group
results = results |> ungroup()
results <- results %>%
  group_by(group, part) %>%
  mutate(R2_norm_env = R2 / sum(R2))

results$PFT = factor(
  results$group,
  levels = c(4,1,2,3,5),
  labels = c("Slow","Fast","Short-lived breeder (SLB)","Long-lived pioneer (LLP)","Intermediate")
)
colors2 <- c(
  "Prec" = "#89bde6",       # Pale Blue
  "SR_kW_m2" = "#ffbb78",   # Pale Orange
  "RH_prc" = "#98df8a",     # Pale Green
  "T_max" = "#ff9896",      # Pale Red
  "T_min" = "#c5b0d5"       # Pale Purple
)

results$Environment = results$part

plt =
  ggplot(results, aes(x = PFT, y = R2, fill = Environment)) + #reorder(group, R2)
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = colors2,
    name = names(colors2)
  )+
  facet_wrap(~env) +
  labs(x = "PFT", y = "Explained Variance [%]") +
  theme_minimal() +
  ggthemes::theme_base()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plt

# save plt as png
ggsave("gfoe2024/BCI_PFT_env_plot.png", plt, width = 34, height = 15, units = "cm")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## continue simulation with fitted parameters ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

pars_sim = pars[length(pars)]

getPars = function(internalPar, parRange) {
  if (is.vector(parRange)) {
    # Case where internalPar is a 1D tensor and parRange is a vector
    Npar <- length(parRange) / 2
    lower <- parRange[1:Npar]
    upper <- parRange[(Npar + 1):(2 * Npar)]

    out <- torch::torch_sigmoid(internalPar) * (upper - lower) + lower
  } else {
    # Case where internalPar is a matrix and parRange is a matrix
    Npar <- ncol(internalPar)
    out <- list()
    for (i in 1:Npar) {
      lower <- parRange[i, 1, drop = FALSE]
      upper <- parRange[i, 2, drop = FALSE]
      out[[i]] <- torch::torch_sigmoid(internalPar[, i, drop = FALSE]) * (upper - lower) + lower
    }
    out <- torch::torch_cat(out, dim = 2L)
  }
  return(out)
}

speciesPars_ranges = list(
  parGrowth = rbind(
    c(0.01, 0.99),
    c(0.01, 4)
  ),
  parMort = rbind(
    c(0.01, 0.99),
    c(0, 4)
  ),
  parReg = c(0.01, 0.99),
  parHeight = c(0.3, 0.7),
  parGrowthEnv = rbind(
    c(-1, 1),
    c(-1, 1)
  ),
  parMortEnv = rbind(
    c(-2, 2),
    c(-2, 2)
  ),
  parRegEnv = rbind(
    c(-2, 2),
    c(-2, 2)
  ))

par_backtr = list()
# i=names(pars_sim[[1]])[5]
for(i in names(pars_sim[[1]])) {
  # check if i beginns with "nn" as characetrs
  if(substr(i, 1, 2) == "nn") {
    par_backtr[[i]] = pars_sim[[1]][[i]]
  }else{
    par_backtr[[i]] = torch::as_array(getPars(internalPar = pars_sim[[1]][[i]], parRange = speciesPars_ranges[[i]]))
  }
}

obs_df = fread("data/calibration-data/BCI-1h-patch-V2/obs_df.csv")
env = fread("data/calibration-data/BCI-1h-patch-V2/env_dt.csv")
df = fread("data/calibration-data/BCI-1h-patch-V2/stand_dt.csv")

missing_years =
  lapply(1982:1984, function(year) {
    tmp = env[2,]
    tmp$year = year
    return(tmp)
  }) |> rbindlist()



env = rbindlist(list(missing_years, env[-1,]))
env = env[!year %in% 2016:2019]

env =
  lapply(unique(df$siteID), function(site) {
    tmp = env
    tmp$siteID = site
    return(tmp)
  }) |> rbindlist()

env$Prec = scale(env$Prec)
env$T_max = scale(env$T_max)
env$T_min = scale(env$T_min)
env$SR_kW_m2 = scale(env$SR_kW_m2)
env$RH_prc = scale(env$RH_prc)
df$AL = as.numeric(df$AL)
str(df)

# env = env[siteID == 1]
# obs_df2 = rbindlist(lapply(1:10, function(i) {
#   obs_df$patchID = i
#   return(obs_df)
# }))

# sp= cohort1$species_r
# sp[sp==0] = 1L
# sp[is.na(sp)] = 1L
# cohort2 = FINN::CohortMat$new(dbh = cohort1$dbh_r, trees = cohort1$trees_r, species = sp, sp = 195)
# pdat = data.table()
# for(i in 1:10) {
iters_start = seq(1, 400, 100)
iters_end = seq(100, 500, 100)
out_dt = data.table()
for(i in length(iters_start)){
  cohort1 <- FINN::CohortMat$new(obs_df[siteID %in% iters_start[i]:iters_end[i]], sp = 195)
  env_tmp = env[siteID %in% iters_start[i]:iters_end[i]]
  orig_siteID = unique(env_tmp$siteID)
  env_tmp[, siteID := as.integer(as.character(factor(siteID, levels = orig_siteID, labels = 1:100))),]
  predictions =
    simulateForest(env = env[siteID %in% iters_start[i]:iters_end[i]],
                   height = as.vector(par_backtr$parHeight),
                   sp = nrow(par_backtr$parReg),
                   init = cohort1,
                   patch_size = 1,
                   patches = 1,
                   growthProcess = createProcess(~., func = growth, initEnv = par_backtr$nnGrowth.0.weight,initSpecies = par_backtr$parGrowth),
                   mortalityProcess = createProcess(~., func = mortality, initEnv = par_backtr$nnMort.0.weight,initSpecies = par_backtr$parMort),
                   regenerationProcess = createProcess(~., func = regeneration, initEnv = par_backtr$nnReg.0.weight,initSpecies = as.vector(par_backtr$parReg)),
                   device = "cpu")

  tmp_dt <- predictions$wide$site
  tmp_dt$rep = i
  env_tmp[, siteID := as.integer(as.character(factor(siteID, levels = orig_siteID, labels = 1:100))),]
  out_dt <- rbind(out_dt, tmp_dt)
}
  # pdat = rbind(pdat, tmp_dt)
# }
pdat[,year := year+1984,]
pdat <- merge(df, pdat, by = c("siteID", "year", "species"))
library(ggplot2)
ggplot() +
  geom_line(aes(x = year, y = ba, group = interaction(rep,species), color = factor(species)), data = pdat[, .(ba = sum(ba.y)), by = .(year, rep, species, siteID)], alpha = 0.3) +
  geom_point(aes(x = year, y = ba, color = factor(species)), data = pdat[, .(ba = sum(ba.x)/100), by = .(year, siteID, species)])+
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  theme(legend.position = "none")+
  facet_wrap(~siteID)







