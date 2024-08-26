# read RDS file
pars <- readRDS("data/calibration-output/parameters2.RDS")

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
kmeans_result1axis <- kmeans(pca_result$x[, 1], centers = 3, nstart = 25)

# Create a data frame with the PCA results and cluster assignments
pca_data <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2],
                       cluster = as.factor(kmeans_result$cluster),
                       cluster1axis = as.factor(kmeans_result1axis$cluster)
                       )

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

assigned_species <- fread("data/calibration-data/BCI/species_assigned.csv")
ruger_species <- data.table(readxl::read_xlsx("data/calibration-output/aaz4797_ruger_data_s1.xlsx", sheet = 2))

pca_data$speciesID = as.integer(rownames(pca_data))

pca_data <- merge(pca_data, assigned_species, by = "speciesID")
ruger_species$sp <- tolower(ruger_species$sp)

pca_data2 <- data.table(merge(pca_data, ruger_species[,.(sp, PFT_1axis, PFT_2axes, PC1score, PC2score)], by = "sp", all = T))

model_dat <- pca_data2[!is.na(PC1score)][order(PC1score)]
plot(PC1 ~ PC1score, model_dat, xlab = "Rüger PC1", ylab = "FINN PC1")
fit <- lm(PC1 ~ PC1score, model_dat)
pred <- predict(fit, newdata = model_dat, interval = "confidence")
lines(model_dat$PC1score, pred[,1], col = "red")
# add confidence interval
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


library(data.table)
library(ggplot2)

# Remove rows with NAs in either cluster or PFT_2axes
pca_data2 <- pca_data2[!is.na(cluster) & !is.na(PFT_2axes)]

# Create a contingency table
contingency_table <- table(pca_data2$cluster, pca_data2$PFT_2axes)

# Find the best matching pairs between the cluster and PFT_2axes
matching <- apply(contingency_table, 1, function(row) {
  colnames(contingency_table)[which.max(row)]
})

# Map the clusters based on the best matching pairs
pca_data2[, aligned_cluster := factor(cluster, labels = matching)]

# Count how many species were assigned to the same group
matching_entries <- pca_data2[aligned_cluster == PFT_2axes, .N]

# Calculate the total number of species
total_entries <- nrow(pca_data2)

# Calculate the percentage of species assigned to the exact same group
percentage_matching <- (matching_entries / total_entries) * 100

# Print the percentage
print(paste("Percentage of species assigned to the exact same group:", round(percentage_matching, 2), "%"))

# Visualization: Create a bar plot to show the comparison
ggplot(pca_data2, aes(x=aligned_cluster, fill=PFT_2axes)) +
  geom_bar(position="dodge") +
  labs(title="Comparison of Classification Methods",
       x="Cluster (Method 1)",
       y="Count",
       fill="Cluster (Method 2)") +
  theme_minimal()

# Alternatively, you can visualize the contingency table as a heatmap
ggplot(as.data.frame(contingency_table), aes(Var1, Var2, fill=Freq)) +
  geom_tile() +
  labs(title="Heatmap of Classification Agreement",
       x="Cluster (Method 1)",
       y="Cluster (Method 2)",
       fill="Number of Species") +
  scale_fill_gradient(low="white", high="red") +
  theme_minimal()

# Create the plot
ggplot(pca_data2, aes(x = PC1, y = PC2, color = factor(PFT_2axes))) +
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




library(data.table)
library(mclust)

# Remove rows with NAs in either cluster or PFT_2axes
pca_data2 <- pca_data2[!is.na(cluster) & !is.na(PFT_2axes)]

# Ensure that the clusters are factors
pca_data2[, cluster := as.factor(cluster)]
pca_data2[, cluster1axis := as.factor(cluster1axis)]
pca_data2[, PFT_2axes := as.factor(PFT_2axes)]
pca_data2[, PFT_1axis := as.factor(PFT_1axis)]

# Calculate the Adjusted Rand Index
ari2 <- adjustedRandIndex(pca_data2$cluster, pca_data2$PFT_2axes)
ari1 <- adjustedRandIndex(pca_data2$cluster1axis, pca_data2$PFT_1axis)

# Print the Adjusted Rand Index
print(paste("Adjusted Rand Index:", round(ari2, 4)))
print(paste("Adjusted Rand Index:", round(ari1, 4)))




