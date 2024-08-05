#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose: create plots for the roxygen documentation of each function
## Date: Mon Aug  5 13:56:54 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=


create_help_plots <- function(){
  dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
  library(ggplot2)

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
}


