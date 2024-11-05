library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env
iterations = 10
all_pred <- data.table()
for(i in 1:iterations){


  Nsp = 10
  shadeSP = runif(Nsp, 0.01, 0.99)
  Npatches = 50
  Nsites = 1
  Ntimesteps = 500

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
  # dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.001)*1
  dist_dt$intensity = 0
  # dist_dt$intensity[250] = 1
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
ggplot(pred$site[, .(value = mean(value)), by = .(year, species, variable)],
       aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year", y = "Value") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_discrete(name = "Species")+
  facet_wrap(~variable, scales = "free_y", ncol = 1, strip.position = "left") +
  theme(axis.title.y = element_blank(),
    strip.placement = "outside",  # Places the facet labels outside the plotting area
    strip.text.y.left = element_text(angle = 90)  # Ensures the facet labels are horizontal
  )

library(gganimate)

# List to store individual animation plots
animations <- list()

# Iterate through each iteration value in your data
for (iter in unique(all_pred$iteration)) {
  # Filter data for the current iteration
  p_dat = all_pred[iteration == iter, .(value = mean(value)), by = .(year, species, variable, iteration)]
  # Assuming `pred$site` has an 'iteration' column
  animated_plot <- ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
    geom_line() +
    labs(x = "Year", y = "Value") +
    theme_minimal() +
    coord_cartesian(ylim = c(0, NA)) +
    facet_wrap(~variable, scales = "free_y", ncol = 1, strip.position = "left") +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.text = element_text(size = 15), legend.text = element_text(size = 15),legend.title = element_text(size = 20),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, size = 20),
    ) +
    # increase linewidth in legend
    guides(color = guide_legend(title = "Species", override.aes = list(linewidth = 5))) +
    scale_color_discrete(name = "Species")+
    transition_states(iteration, transition_length = 3, state_length = 1) +  # Transition between iterations
    transition_reveal(year) +  # Draw the line over time within each iteration
    ease_aes('cubic-in-out')  # Smooth transition

  gif_filename <- paste0("iteration_", iter, ".gif")
  anim_save(gif_filename, animate(animated_plot, nframes = 100, fps = 10, width = 500, height = 1000))

  # Store the animation in the list
  animations[[iter]] <- gif_filename
  # Render and save the plot as a GIF
  # anim_save("playground/animated_plot.gif", animate(animated_plot, nframes = 20, fps = 5, width = 500, height = 1000))
}

library(magick)

# Read all the GIFs
gif_images <- image_read(unlist(animations))

# Combine them into a single animation
combined_gif <- image_join(gif_images) %>%
  image_animate(fps = 10)  # Control the frame rate for the combined GIF

# Save the combined GIF
image_write(combined_gif, "combined_iterations.gif")
