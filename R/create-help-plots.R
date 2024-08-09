#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose: create plots for the roxygen documentation of each function
## Date: Mon Aug  5 13:56:54 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=


create_help_plots <- function(){
  dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
  library(ggplot2)
  library(viridis)
  library(data.table)

  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  ## height ####
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  test_dt =
    data.frame(
      expand.grid(
        list(
          dbh = seq(0, 500,1),
          parHeight = seq(0,1,0.1)
        )
      )
    )
  test_dt$rowID = 1:nrow(test_dt)

  test_dt$height = height(test_dt$dbh, test_dt$parHeight)

  p_height <-
    ggplot(test_dt, aes(x = dbh, y = height, group = parHeight,color = parHeight))+
    ylab("height(dbh,parHeight)")+
    xlab("dbh [cm]")+
    geom_line()+
    scale_y_continuous(breaks = seq(0, 200, 10))+
    geom_label(data = test_dt[which(test_dt$dbh == max(test_dt$dbh)),], aes(label = parHeight, color = parHeight))

  ggsave("man/figures/height_plot1.png", p_height, width = 13, height = 10, units = "cm")

  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  ## BA_stand ####
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  # Define vectors
  trees_vec <- c(0:500, 10^(seq(2, 4, length.out = 20)))
  dbh_vec <- seq(1, 200, 1)
  patch_size_vec <- c(0.01, 0.05, 0.08, seq(0.1, 1, 0.2), 1)

  # Initialize an empty data frame
  test_dt_out <- data.frame()

  # Loop through patch sizes and generate test data
  for (i in patch_size_vec) {
    cohort_df1 <- expand.grid(
      trees_ha = trees_vec,
      patch_size_ha = i,
      dbh = dbh_vec
    )

    cohort_df1 <- data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      cohort_df1
    )

    cohort_df1$siteID <- 1:nrow(cohort_df1)
    cohort_df1$trees <- round(cohort_df1$trees_ha * i)

    cohort <- CohortMat$new(obs_df = cohort_df1)

    basal_area <- BA_stand(cohort$dbh, cohort$trees, patch_size_ha = i)
    light <- competition(dbh = cohort$dbh, species = cohort$species,
                         trees = cohort$trees, parHeight = torch::torch_tensor(c(0.5)),
                         patch_size_ha = i,
                         h = 0)

    cohort_df1$basal_area <- torch::as_array(basal_area)[, 1, 1]
    cohort_df1$light <- torch::as_array(light)[, 1, 1]

    test_dt_out <- rbind(test_dt_out, cohort_df1)
  }


  cohort_df1_plot <- test_dt_out[test_dt_out$patch_size == 0.1 & test_dt_out$basal_area > 0,]
  p_BA_stand2 <- ggplot(cohort_df1_plot, aes(x = dbh, y = basal_area, color = trees, group = trees)) +
    geom_line() +
    ylab("Basal Area (m^2/ha)") +
    xlab("Diameter at Breast Height (cm)") +
    coord_cartesian(ylim = c(0,59.999))+
    scale_y_continuous(breaks = seq(0, 60, 10)) +
    # add a multicolor left log10 scaled gradient for the color scale
    scale_color_viridis_c(name = "Trees per ha", trans = "log10", option = "magma", direction = -1)

  # Cap the basal_area at max_ba
  max_ba <- 60
  test_dt_out$basal_area[test_dt_out$basal_area > max_ba] <- max_ba

  p_BA_stand1 <- ggplot(test_dt_out[test_dt_out$trees_ha <= 500, ], aes(y = factor(trees_ha), x = dbh, fill = basal_area)) +
    geom_tile() +
    ylab("trees [N/ha]") +
    xlab("mean stand dbh [cm]") +
    guides(fill = guide_legend(title = "basal area\n[m^2/ha]")) +
    # ggtitle("BA_stand", paste0(deparse(BA_stand), collapse = "\n")) +
    scale_fill_viridis_c(breaks = seq(0, max_ba, 10), labels = c(seq(0, max_ba - 10, 10), paste0(">", max_ba))) +
    scale_y_discrete(breaks = seq(0, 500, 100)) +
    facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ", x, " ha")))

  ggsave("man/figures/BA_stand_plot1.png", p_BA_stand1, width = 15, height = 10, units = "cm")
  ggsave("man/figures/BA_stand_plot2.png", p_BA_stand2, width = 15, height = 10, units = "cm")

  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  ## competition ####
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  # Generate plots for inclusion in the documentation

  # Define vectors
  trees_vec <- c(0:1000)
  # dbh_vec <- seq(1, 10, 1)
  dbh_vec <- seq(0, 500, 1)
  patch_size_ha <- 1

  # Sensitivity of competition to parHeight (species)
  parHeight_vec <- seq(0, 1, 0.1)
  cohort_dt_out <- data.frame()
    i=1
  for (i in 1:length(parHeight_vec)) {
    parHeight_i = parHeight_vec[i]
    cohort_df1 <- data.frame(
      patchID = 1,
      cohortID = i,
      species = i,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = 0.1,
        dbh = dbh_vec
      )
    )
    cohort_df1$siteID <- 1:nrow(cohort_df1)
    cohort_df1$parHeight <- parHeight_i
    cohort_dt_out <- rbind(cohort_dt_out, cohort_df1)
  }


    cohort_dt_out$trees <- round(cohort_dt_out$trees_ha * patch_size_ha)
    cohort <- CohortMat$new(obs_df = cohort_dt_out)

    basal_area <- BA_stand(cohort$dbh, cohort$trees, patch_size_ha = patch_size_ha)
    # cohortHeights2 = height(cohort$dbh, torch::torch_tensor(parHeight_vec))$unsqueeze(4)
    # dim(torch::as_array(cohortHeights2))
    # cohort_dt_out$height <- c(torch::as_array(cohortHeights2)[,,1,], torch::as_array(cohortHeights2)[,,2,])
    light <- competition(dbh = cohort$dbh, species = cohort$species,
                         trees = cohort$trees, parHeight = torch::torch_tensor(parHeight_vec),
                         h = NULL, patch_size_ha = patch_size_ha)
    cohort_dt_out$basal_area <- as.vector(torch::as_array(basal_area)[, 1, ])
    cohort_dt_out$light <- as.vector(torch::as_array(light)[, 1, ])

    # average light over round(basal_area), species, parHeight with base r
    cohort_dt_out <- data.table(cohort_dt_out)
    cohort_dt_out <- cohort_dt_out[, .(light = mean(light)), by = .(species, parHeight, basal_area = round(basal_area))]


    p_competition1 <- ggplot(
      cohort_dt_out[cohort_dt_out$basal_area <= 60, ],
      aes(x = basal_area, y = light, color = parHeight, group = interaction(species,parHeight))) +
      # geom_tile() +
      geom_line()+
      # coord_cartesian(ylim = c(0.8,1), xlim = c(0,10))+
      # ylab("trees [N/ha]") +
      # ylab("light [%]") +
      # xlab("mean stand dbh [cm]") +
      scale_color_viridis_c(name = "species' parHeight",
                            option = "cividis",
                            direction = 1,
                            breaks = seq(0,1,0.1),
                            labels = c(
                              "0.0 = forest floor",
                              "0.1 = small tree species",
                              seq(0.2,.9,0.1),
                              "1.0 = tall tree species"
                              )) +
    xlab("total stand basal area [m^2/ha]")+
    theme_classic()+
    theme(
      legend.position = c(0.8, 0.6),  # Adjust coordinates as needed
      legend.background = element_rect(fill = NA)
    )
    p_competition1
    ggsave(paste0("man/figures/competition_parHeight.png"), p_competition1, width = 15, height = 10, units = "cm")

  # Sensitivity of competition to h: light availability vs. height
  patch_size <- 1
  trees_vec <- c(1:500)
  dbh_vec <- seq(1, 200, 10)

  cohort_df1 <- data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees = round(trees_vec * patch_size),
      dbh = dbh_vec
    )
  )
  cohort_df1$siteID <- 1:nrow(cohort_df1)
  cohort <- CohortMat$new(obs_df = cohort_df1)

  light <- competition(cohort$dbh, cohort$species, cohort$trees,
                       parHeight = torch::torch_tensor(0.5), h = 0, patch_size_ha = patch_size)
  cohort_df1$light <- torch::as_array(light)[, 1, 1]

  # use BA_stand to calculate basal area of the stand
  basal_area <- BA_stand(cohort$dbh, cohort$trees, patch_size_ha = patch_size)
  cohort_df1$basal_area <- torch::as_array(basal_area)[, 1, 1]

  # cap basal area at 60 m^2/ha
  cohort_df1$basal_area[cohort_df1$basal_area > 60] <- 60
  # cohort_df1[cohort_df1$basal_area > 60]

  p_competition2 <-    ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light*100))+
    geom_tile()+
    ylab("trees/ha [N]")+
    xlab("mean stand dbh [cm]")+
    guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
    ggtitle("competition",paste0(deparse(competition),collapse = "\n"))
  p_competition2

  ggsave("man/figures/competition_light.png", p_competition2, width = 15, height = 10, units = "cm")


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




}


