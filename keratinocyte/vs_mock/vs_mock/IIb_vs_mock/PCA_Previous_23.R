# Load required libraries (if not already loaded)
# install.packages("ggplot2")  # Uncomment and run if you haven't installed ggplot2
library(ggplot2)

# Assuming your count matrix is stored in a variable called 'count_matrix'
# Replace 'count_matrix' with your actual variable name
Xx<-read.csv("vsd_raw.csv")
# Extract the columns containing sample data
sample_data <- Xx[, 2:13]

# Standardize the data (mean center and scale to unit variance)
sample_data_scaled <- scale(sample_data)

# Perform PCA
pca_result <- prcomp(sample_data_scaled)

# Plot the variance explained by each principal component
summary(pca_result)

# Create a scree plot to visualize the explained variance
scree_data <- data.frame(
  PC = 1:ncol(sample_data_scaled),
  VarianceExplained = (pca_result$rotation^2) / sum(pca_result$rotation^2)
)

# Plot the scree plot
ggplot(data = scree_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity") +
  labs(title = "Scree Plot", x = "Principal Component (PC)", y = "Variance Explained")

# You can also visualize the PCA results with a biplot
biplot(pca_result)

# Access the scores and loadings of the principal components
scores <- pca_result$x  # Scores (sample projections onto PCs)
loadings <- pca_result$rotation  # Loadings (correlation between variables and PCs)





# Assuming you already have your 'sample_data_scaled' from the previous code

# Perform PCA
pca_result <- prcomp(sample_data_scaled)

# Create a PCA plot with only sample data as dots
plot(pca_result$rotation[,1], pca_result$rotation[,2], pch = 19, col = "blue", xlab = "PC1", ylab = "PC2", main = "PCA Plot")



# Assuming you already have your 'sample_data_scaled' from the previous code

# Perform PCA
pca_result <- prcomp(sample_data_scaled)

# Create a PCA plot with labeled sample data points
plot(pca_result$rotation[, 1], pca_result$rotation[, 2], pch = 19, col = "blue", xlab = "PC1", ylab = "PC2", main = "PCA Plot")

# Add text labels to the data points
text(pca_result$rotation[, 1], pca_result$rotation[, 2], labels = colnames(sample_data), pos = 3, cex = 0.7)

# Assuming you already have your 'sample_data_scaled' from the previous code

# Create a vector of colors for the four groups
group_colors <- c("red", "blue", "green", "purple")

# Perform PCA
pca_result <- prcomp(sample_data_scaled)

# Create a PCA plot with labeled sample data points, colored by groups
plot(pca_result$rotation[, 1], pca_result$rotation[, 2], pch = 19, col = group_colors, xlab = "PC1", ylab = "PC2", main = "PCA Plot of Keratinocyte Raw Count")

# Add text labels to the data points
text(pca_result$rotation[, 1], pca_result$rotation[, 2], labels = colnames(sample_data), pos = 3, cex = 0.7)
# Create a legend
legend("topleft", legend = c("Mock", "IIa", "IIb", "I"), fill = group_colors)

