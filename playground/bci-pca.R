# read RDS file
pars <- readRDS("data/calibration-output/parameters.RDS")

length(pars)

pars2 = pars[[length(pars)]]

out <- data.table()
i=2
for(i in 1:length(pars2)){
  par_name = names(pars2)[i]
  print(i)
  print(par_name)
  dt <- data.table(pars2[[i]])
  colnames(dt) <- paste0(par_name, "_", 1:ncol(pars2[[i]]))
  out = cbind(out, dt)
}
pars2$parMort

# only select columns with "par" in the name
data = out[,grep("par", colnames(out)), with = FALSE]
# data = out


# Load necessary libraries
library(ggplot2)

# Assuming 'data' is already loaded as shown in your screenshot
# If 'data' is not yet a data frame, you need to load or prepare it accordingly

# Perform PCA on the data
pca_result <- prcomp(data, scale. = TRUE)

# Perform K-means clustering on the first two principal components
set.seed(123)
kmeans_result <- kmeans(pca_result$x[, 1:2], centers = 5, nstart = 25)

# Create a data frame with the PCA results and cluster assignments
pca_data <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2],
                       cluster = as.factor(kmeans_result$cluster))

# Extract loadings for the first two principal components
loadings <- pca_result$rotation[, 1:2]
loadings_df <- data.frame(Variable = rownames(loadings),
                          PC1 = loadings[, 1] * 5,  # Scaling for better visualization
                          PC2 = loadings[, 2] * 5)

# Create the plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("magenta", "cyan", "green", "blue", "orange"),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3",
                                "Cluster 4", "Cluster 5")) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = Variable),
            hjust = 0.5, vjust = 0.5, color = "black")
