library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)


#set the working directory
setwd("D:/DEG/GSE219036_6_june/keratinocyte")

#load the count data
count_data <- read.csv("count_keratinocyte.csv", header=TRUE,row.names = 1)
colnames(count_data)
head(count_data)

#load the sample info
sample_info <- read.csv("meta_keratinocyte.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Cell_type <- factor(sample_info$Cell_type)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)


unique(sample_info$Treatment)
unique(sample_info$Cell_type)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ Treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)

#set the referene for the treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c ("mock","MPXV clade I infected", "MPXV clade IIb infected", "MPXV clade IIa infected")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the factor level
dds$Treatment <- relevel(dds$Treatment , ref = "mock")
dds$Treatment


#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)



#for MPXV clade IIa infected
deseq_results_IIa <- results(dds, contrast=c("Treatment", "MPXV clade IIa infected", "mock"))
deseq_results_IIa <- as.data.frame(deseq_results_IIa)
head(deseq_results_IIa)
deseq_results_IIa_0.05 <- results(dds,contrast=c("Treatment", "MPXV clade IIa infected", "mock"), alpha = 0.05)
deseq_results_IIa_0.05 <- as.data.frame(deseq_results_IIa_0.05)
head(deseq_results_IIa_0.05)

# create histogram plot of p-values
hist(deseq_results_IIa_0.05$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE219036 Frequencies of padj-values")

#for MPXV clade IIb infected
deseq_results_IIb <- results(dds, contrast=c("Treatment", "MPXV clade IIb infected", "mock"))
deseq_results_IIb <- as.data.frame(deseq_results_IIb)
head(deseq_results_IIb)
deseq_results_IIb_0.05 <- results(dds,contrast=c("Treatment", "MPXV clade IIb infected", "mock"), alpha = 0.05)
deseq_results_IIb_0.05 <- as.data.frame(deseq_results_IIb_0.05)
head(deseq_results_IIb_0.05)

# create histogram plot of p-values
hist(deseq_results_IIb_0.05$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE219036 Frequencies of padj-values")

