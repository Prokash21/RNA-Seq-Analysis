## R package installer for this project
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
cran_pkgs <- c(
  "tidyverse",
  "pheatmap",
  "ComplexHeatmap",
  "RColorBrewer",
  "cowplot",
  "EnhancedVolcano"
)
install.packages(setdiff(cran_pkgs, rownames(installed.packages())))

# Bioconductor packages
bioc_pkgs <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "tximport",
  "GenomicFeatures",
  "SummarizedExperiment",
  "BiocParallel"
)
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), ask = FALSE)

message("R package installation completed. Restart R or RStudio if packages were updated.")
