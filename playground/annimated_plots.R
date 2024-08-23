library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env
iterations = 5
Ntimesteps = 1000
dist_intensities = c(0.2,0.4,0.6)
dist_years = c(250,400,700)
all_pred <- data.table()
FINN.seed(1)
for(i in 1:iterations){


  Nsp = 10
  Npatches = 50
  Nsites = 1
  if(i != 2){
    shadeSP = runif(Nsp, 0.01, 0.99)
    speciesPars = list(
      speciesID = 1:Nsp,
      parGrowth = matrix(c(
        shadeSP,
        runif(Nsp, 1, 2)
        ),Nsp, 2),
      parMort = matrix(c(
        shadeSP,
        runif(Nsp, 1, 4)
        ),Nsp, 2),
      parReg = shadeSP,
      parHeight = runif(Nsp, 0.3, 0.7),
      parGrowthEnv = list(matrix(c(
        runif(Nsp, -1, 1),
        runif(Nsp, -1, 1)
        ),Nsp, 2)),
      parMortEnv = list(matrix(c(
        runif(Nsp, 0, 0.5),
        runif(Nsp, 0, .5)
        ), Nsp, 2)),
      parRegEnv = list(matrix(c(
        runif(Nsp, 0, 2),
        runif(Nsp, -2, 2)
        ),Nsp, 2))
    )
  }
  site_dt <- data.table(
    expand.grid(
      list(
        siteID = 1:Nsites,
        year = 1:Ntimesteps
      )
    )
  )

  # hist(rbeta(Ntimesteps*Nsites, 1, 10), xlim = c(0,1))
  dist_dt <- site_dt
  # dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.1)*
  dist_dt$intensity = 0
  if(i > 1){
    dist_dt[year == dist_years[1], intensity := dist_intensities[1]]
    dist_dt[year == dist_years[2], intensity := dist_intensities[2]]
    dist_dt[year == dist_years[3], intensity := dist_intensities[3]]
  }
  # dist_dt$intensity[500] = 0.5
  # dist_dt$intensity[750] = 0.25
  # dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.01)*rbeta(Ntimesteps*Nsites, 1, 5)

  env_dt <- site_dt
  env_dt$env1 = rep(0, Ntimesteps)
  # env_dt$env1 = rep(seq(-2,2,length.out = Nsites), Ntimesteps)

  system.time({
    predictions =
      simulateForest(env = env_dt,
                     disturbance = dist_dt,
                     sp = Nsp,
                     patches=Npatches,
                     # height = speciesPars$parHeight,speciesPars_ranges = speciesPars_ranges,
                     # growthProcess = createProcess(~1+env1, func = growth, initEnv = list(matrix(c(runif(),3,4), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                     # mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = list(matrix(c(1, -3), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                     # regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = list(matrix(c(4, 0), 1, 2)), initSpecies = matrix(c(0.2), Nsp)),
                     growthProcess = createProcess(~1+env1, func = growth, initEnv = speciesPars$parGrowthEnv, initSpecies = speciesPars$parGrowth),
                     mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = speciesPars$parMortEnv, initSpecies = speciesPars$parMort),
                     regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = speciesPars$parRegEnv, initSpecies = speciesPars$parReg),
                     device = "cpu", debug = F)
  })

  pred = predictions$long$site
  pred$iteration = i
  all_pred = rbind(all_pred, pred)
  cat(i, "of",iterations, "\n")
}

library(ggplot2)
# ggplot(all_pred[iteration == 1, .(value = mean(value)), by = .(year, species, variable)],
#        aes(x = year, y = value, color = factor(species))) +
#   geom_line() +
#   labs(x = "Year", y = "Value") +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, NA)) +
#   scale_color_discrete(name = "Species")+
#   facet_wrap(~variable, scales = "free_y", ncol = 1, strip.position = "left") +
#   theme(axis.title.y = element_blank(),
#     strip.placement = "outside",  # Places the facet labels outside the plotting area
#     strip.text.y.left = element_text(angle = 90)  # Ensures the facet labels are horizontal
#   )

library(gganimate)

if(!dir.exists("vignettes/annimated-plots")){
  dir.create("vignettes/annimated-plots", recursive = T)
}

all_pred[,value2 := value,]

# calculate ba, reg, and trees per ha
all_pred[variable %in% c("ba", "trees", "reg"), value2 := value/0.1,]
# calculate AL and mort as %
all_pred[variable %in% c("AL", "mort"), value2 := round(value*100,4),]

unique(all_pred$variable)
# rename variable names with units and more descriptive names
# dbh    ba     trees  AL     growth mort   reg
# cm     m^2/ha trees/ha m/yr   1/yr   1/yr  1/yr
# Update variable2 labels with proper mathematical notation
all_pred[variable == "dbh", variable2 := "DBH~(cm)"]
all_pred[variable == "ba", variable2 := "Basal~Area~(m^2/ha)"]
all_pred[variable == "trees", variable2 := "N~of~Trees~(1/ha)"]
all_pred[variable == "AL", variable2 := "Available~Light~('%')"]
all_pred[variable == "growth", variable2 := "DBH~Growth~(cm)"]
all_pred[variable == "mort", variable2 := "Mortality~('%')"]
all_pred[variable == "reg", variable2 := "Regeneration~(1/ha)"]

# Set variable2 as a factor and order it
all_pred[, variable2 := factor(
  variable2,
  levels = c(
    "Basal~Area~(m^2/ha)",
    "DBH~Growth~(cm)",
    "N~of~Trees~(1/ha)",
    "Mortality~('%')",
    "DBH~(cm)",
    "Regeneration~(1/ha)",
    "Available~Light~('%')"
  ), ordered = TRUE
)]

# Define the shift_legend3 function (no changes needed here)
shift_legend3 <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) |> gtable::gtable_filter("panel") |>
    with(setNames(grobs, layout$name)) |> purrr::keep(~identical(.x,zeroGrob()))

  if (length(pnls) == 0) stop("No empty facets in the plot")

  lemon::reposition_legend(p, "center", panel=names(pnls))
}

# Prepare the data for plotting
p_dat <- all_pred[iteration == 1, .(value = mean(value2)), by = .(year, species, variable2, iteration)]

# Create the ggplot
p <- ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year", y = "Value") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left", labeller = label_parsed) +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  guides(color = guide_legend(title = "Species", override.aes = list(linewidth = 5), ncol = 3, title.position = "left")) +
  scale_color_discrete(name = "Species")

# To potentially shift the legend, use shift_legend3 function
# p <- shift_legend3(p)

# Print the plot
# print(p)

# shift_legend3(p)
# List to store individual animation plots
animations <- list()

# Iterate through each iteration value in your data
iter = 1
for (iter in unique(all_pred$iteration)) {
  # Filter data for the current iteration
  p_dat = all_pred[iteration == iter, .(value = mean(value2)), by = .(year, species, variable = variable2, iteration)]
  # Assuming `pred$site` has an 'iteration' column
  if(iter == 1) years_plot = 0 else years_plot = 1
  animated_plot <- ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
    geom_line() +
    labs(x = "Year", y = "Value") +
    theme_minimal() +
    coord_cartesian(ylim = c(0, NA)) +
    facet_wrap(~variable, scales = "free_y", ncol = 2, strip.position = "left", labeller = label_parsed) +
    theme(
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90)
    ) +
    geom_vline(xintercept = dist_years, linetype = "dashed", linewidth = 0.4, alpha = years_plot)+
    guides(color = guide_legend(title = "Species", override.aes = list(linewidth = 5), ncol = 3, title.position = "left")) +
    scale_color_discrete(name = "Species")+
    transition_states(iteration, transition_length = 3, state_length = 1) +  # Transition between iterations
    transition_reveal(year) +  # Draw the line over time within each iteration
    ease_aes('cubic-in-out')  # Smooth transition

  gif_filename <- paste0("vignettes/annimated-plots/iteration_", iter, ".gif")
  anim_save(gif_filename, animate(animated_plot, nframes = 100, fps = 10, width = 1800, height = 1200, res = 190))

  # Store the animation in the list
  animations[[iter]] <- gif_filename
  # Render and save the plot as a GIF
  # anim_save("playground/animated_plot.gif", animate(animated_plot, nframes = 20, fps = 5, width = 500, height = 1000))
}

library(magick)

# Read all the GIFs
gif_images <- image_read(unlist(animations))

# Combine them into a single animation
combined_gif <- image_join(gif_images) #%>%
  #image_animate(fps = 10)  # Control the frame rate for the combined GIF

# Save the combined GIF
image_write(combined_gif, "vignettes/annimated-plots/combined_iterations.gif")

# Read the GIF
gif <- image_read("vignettes/annimated-plots/combined_iterations.gif")

# # Compress the GIF by reducing colors, resizing, and optimizing
# gif <- gif %>%
#   image_quantize(max = 128)  # Reduce number of colors
#   # image_scale("50%") %>%         # Resize to 50% of the original dimensions
#   # image_optimize()               # Optimize the GIF
# # Save the compressed GIF
# image_write(gif, "vignettes/annimated-plots/combined_iterations_compressed1.gif")

gif <- gif %>%
  # image_quantize(max = 128)  # Reduce number of colors
  image_scale("50%")         # Resize to 50% of the original dimensions
  # image_optimize()               # Optimize the GIF
# Save the compressed GIF
image_write(gif, "vignettes/annimated-plots/combined_iterations_compressed2.gif")
