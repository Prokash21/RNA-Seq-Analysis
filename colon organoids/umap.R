# Load required packages
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(umap)

# Set the working directory
setwd("D:/DEG/GSE219036_6_june/colon organoids")

# Load count data
count_data <- read.csv("count_colon.csv", header = TRUE, row.names = 1)

# Load sample information
sample_info <- read.csv("meta_colon.csv", header = TRUE, row.names = 1)

# Convert non-integer values to integers in count data
count_data <- round(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~Treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Set the reference for the treatment factor
dds$Treatment <- relevel(dds$Treatment, ref = "MPXV clade I infected")

# Filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Perform statistical tests to identify differentially expressed genes
dds <- DESeq(dds)

# Results and filtering
deseq_results <- results(dds)
deseq_results <- as.data.frame(deseq_results)
filtered <- deseq_results %>% filter(deseq_results$padj < 0.05, abs(deseq_results$log2FoldChange) > 1)

# Normalize counts
normalize_counts <- counts(dds, normalized = TRUE)
############################################################
# UMAP transformation
library(umap)

# Normalize counts (if not already done)
normalize_counts <- counts(dds, normalized = TRUE)

# Check the number of items (rows) in your dataset
num_items <- nrow(normalize_counts)

# Specify an appropriate number of neighbors (e.g., a smaller value)
n_neighbors <- min(10, num_items)  # Adjust the number of neighbors as needed

# UMAP transformation with the specified number of neighbors
umap_result <- umap(t(normalize_counts), n_neighbors = n_neighbors)

# Convert the UMAP result to a data frame
umap_data <- as.data.frame(umap_result$layout)

# Add sample information to the UMAP data frame
umap_data <- cbind(umap_data, sample_info)

# Visualize UMAP results
library(ggplot2)

# Create a scatter plot of UMAP results colored by Treatment
ggplot(umap_data, aes(x = umap_result$layout[, 1], y = umap_result$layout[, 2], color = Treatment)) +
  geom_point() +
  labs(title = "UMAP Visualization") +
  scale_color_brewer(palette = "Set1")  # You can change the palette as needed
#######################################################
