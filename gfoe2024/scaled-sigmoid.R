# Load necessary library
library(viridis)

# Define the scaled sigmoid function
ffun_scaled <- function(environment, parameter, base_steepness = 10) {
  scaled_steepness <- base_steepness / (0.5 - abs(parameter - 0.5))
  ((1 / (1 + exp(-scaled_steepness * (environment - parameter))) - 1 / (1 + exp(scaled_steepness * parameter))) /
         (1 / (1 + exp(-scaled_steepness * (1 - parameter))) - 1 / (1 + exp(scaled_steepness * parameter))))
}

# Define the standard sigmoid function
ffun_sigmoid <- function(environment, parameter, steepness = 10) {
  1 / (1 + exp(-steepness * (environment - parameter)))
}

# Define the step function
ffun_step <- function(environment, parameter) {
  ifelse(environment < parameter, 0, 1)
}

# Function to create and save plots
save_folder = "gfoe2024"
# Create the folder if it doesn't exist
if (!dir.exists(save_folder)) {
  dir.create(save_folder)
}

# Environment levels for plotting
environment <- seq(0, 1, 0.01)

# Parameter values for demonstration
parameters <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# Generate a modern color palette for different parameter values
colors <- viridis(length(parameters))

# Plot and save the scaled sigmoid function
png(file.path(save_folder, "scaled_sigmoid.png"), width = 500, height = 500)
plot(environment, ffun_scaled(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Scaled Sigmoid Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_scaled(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))
dev.off()

# Plot and save the standard sigmoid function
png(file.path(save_folder, "sigmoid.png"), width = 500, height = 500)
plot(environment, ffun_sigmoid(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Standard Sigmoid Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_sigmoid(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))
dev.off()

# Plot and save the step function
png(file.path(save_folder, "step_function.png"), width = 500, height = 500)
plot(environment, ffun_step(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Step Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_step(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))
dev.off()

# Additionally, display the plots in R
par(mfrow = c(1, 3))  # Arrange the plots in a single row

# Display scaled sigmoid function
plot(environment, ffun_scaled(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Scaled Sigmoid Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_scaled(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))

# Display standard sigmoid function
plot(environment, ffun_sigmoid(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Standard Sigmoid Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_sigmoid(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))

# Display step function
plot(environment, ffun_step(environment, parameters[1]), type = "l", ylim = c(0,1),
     col = colors[1], lwd = 2, xlab = "Environment", ylab = "Output", main = "Step Function")
for (i in 2:length(parameters)) {
  lines(environment, ffun_step(environment, parameters[i]), col = colors[i], lwd = 2)
}
legend("right", legend = paste("Parameter =", parameters), col = colors, lty = 1, lwd = 2, bg = rgb(1, 1, 1, 0.3))

# Optionally, specify a different folder
# create_and_save_plots("custom_folder_name")
