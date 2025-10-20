## Bulk RNA-seq R package installer
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

message("Bulk RNA-seq R package installation completed.")
