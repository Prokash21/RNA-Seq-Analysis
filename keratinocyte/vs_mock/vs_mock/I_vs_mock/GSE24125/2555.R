library(affy)
library(GEOquery)
library(tidyverse)
setwd('D:/DEG/GSE219036_6_june/25')
# get supplementary files
getGEOSuppFiles("GSE24125")

# untar files
untar("GSE24125/GSE24125_RAW.tar", exdir = 'data/')

# Load the necessary package for handling compressed files
install.packages("R.utils")
library(R.utils)

# Replace 'GPL10913_Print_1110.adf.gz' with the actual file path
adf_file <- "GPL10913_Print_1110.adf.gz"

# Decompress the ADF file
gunzip(adf_file, remove = TRUE)

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")
ReadAffy("GPL10913_Print_1110.adf.gz")
# performing RMA normalization
normalized.data <- rma(raw.data)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')




library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE148537")

# untar files
untar("GSE148537/GSE148537_RAW.tar", exdir = 'data/')

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")

# performing RMA normalization
normalized.data <- rma(raw.data)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')