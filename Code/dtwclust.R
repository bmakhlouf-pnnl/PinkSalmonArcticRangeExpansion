#install.packages("dtwclust")
library(dtwclust)

# load in the matrix 
dat <- read.csv("/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Data/Cleaned/OtoDataInterpolated_GAM_Matrix.csv")

# Save the first column (labels) before removing it
labels <- sub("\\..*", "", dat[, 1])

# Remove the first column 
dat <- dat[, -1]

################################################################################
pc <- tsclust(dat, type = "partitional", k = 7L, 
              distance = "dtw_basic", centroid = "pam", 
              seed = 3247L, trace = TRUE,
              args = tsclust_args(dist = list(window.size = 20L)))
plot(pc)

################################################################################
hc <- tsclust(dat, type = "hierarchical", k = 10L, 
              distance = "sbd", trace = TRUE,
              control = hierarchical_control(method = "average"))
plot(hc)

################################################################################
##### Run a PCA with the matrix, colored by the label

pca_result <- prcomp(dat, scale. = TRUE)

# Create color mapping for locations
location_colors <- c("Elson" = "#E41A1C", "HART" = "#377EB8", "JACK" = "#4DAF4A")
colors <- location_colors[labels]

# Set up a 1x3 grid for three plots
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

# Plot 1: PC1 vs PC2
plot(pca_result$x[, 1:2], col = colors, 
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
     main = "PC1 vs PC2",
     pch = 19, cex = 3, cex.lab = 1.2, cex.main = 1.4)
legend("topright", legend = names(location_colors), col = location_colors, pch = 19, cex = 1)

# Plot 2: PC1 vs PC3
plot(pca_result$x[, c(1, 3)], col = colors,
     xlab = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
     ylab = paste0("PC3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)"),
     main = "PC1 vs PC3",
     pch = 19, cex = 3, cex.lab = 1.2, cex.main = 1.4)
legend("topright", legend = names(location_colors), col = location_colors, pch = 19, cex = 1)

# Plot 3: PC2 vs PC3
plot(pca_result$x[, c(2, 3)], col = colors,
     xlab = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
     ylab = paste0("PC3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)"),
     main = "PC2 vs PC3",
     pch = 19, cex = 3, cex.lab = 1.2, cex.main = 1.4)
legend("topright", legend = names(location_colors), col = location_colors, pch = 19, cex = 1)

par(mfrow = c(1, 1))

################################################################################
##### 3D PCA Plot with Plotly

library(plotly)

# Create 3D scatter plot
fig <- plot_ly(x = pca_result$x[, 1], 
               y = pca_result$x[, 2], 
               z = pca_result$x[, 3],
               mode = 'markers',
               type = 'scatter3d',
               marker = list(size = 8, color = colors, opacity = 0.8),
               text = labels,
               hoverinfo = 'text') %>%
  layout(
    title = "3D PCA of Otolith Timeseries",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)"))
    )
  )

fig

################################################################################
##### k-NN Classification

library(class)

# Perform k-NN classification with k=5 (using leave-one-out cross-validation)
k <- 5
knn_pred <- knn(dat, dat, labels, k = k)

# Calculate accuracy
accuracy <- sum(knn_pred == labels) / length(labels)
cat("k-NN Classification (k =", k, ")\n")
cat("Accuracy:", round(accuracy * 100, 1), "%\n")
cat("Correct predictions:", sum(knn_pred == labels), "out of", length(labels), "\n\n")

# Confusion matrix
cat("Confusion Matrix:\n")
print(table(Predicted = knn_pred, Actual = labels))