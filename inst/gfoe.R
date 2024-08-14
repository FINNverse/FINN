#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## 2. simulate data ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
sp = 9L
patches = 50L
sites = 1L
initCohort = CohortMat$new(dims = c(sites, patches, 10),
                           dbh = array(1, dim = c(sites, patches, 10)),
                           trees = array(1, dim = c(sites, patches, 10)),
                           sp = sp)




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


env = torch::torch_randn(size = c(sites, 3000L, 2))*0+1
dim(env)
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

melt_dt <- melt(inventory_dt, id.vars = c("siteID","year", "species"), measure.vars = out_names)

ggplot(melt_dt[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = out_names[i]) +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")


# # make animated ggplot that starts at year == 100 and ends at year == 2000
# # load packages for animated plot
# library(gganimate)
# library(magick)
#
# # Create the animated ggplot
# # Create the animated ggplot
# animated_plot <- ggplot(melt_dt[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species), group = species)) +
#   geom_line() +
#   labs(x = "Year",
#        y = "Value") +  # Add dynamic title to show current year
#   theme_minimal() +
#   facet_wrap(~variable, scales = "free") +
#   transition_reveal(year) +               # Animate over the 'year' variable
#   ease_aes('linear')                    # Use linear easing for smooth animation
#
# # Animate and render directly in the plot window and save as a GIF
# anim <- animate(animated_plot, renderer = gifski_renderer(), width = 1000, height = 600, fps = 20)
#
# # Save the animation as a GIF
# anim_save("mtcars_animation.gif", animation = anim)
