library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)


#set the working directory
setwd("D:/DEG/GSE219036_6_june/keratinocyte/vs I/IIa")

#load the count data

count_data <- read.csv("count_keratinocyte.csv", header=TRUE,row.names = 1)
colnames(count_data)
head(count_data)

#load the sample info
sample_info1 <- read.csv("meta_keratinocyte.csv")
sample_info2 <- read.csv("meta_keratinocyte.csv", header = TRUE,row.names = 1)
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
head(dds)

#set the factor level
dds$Treatment <- factor(dds$Treatment, levels = c ("mock","MPXV clade I infected", "MPXV clade IIb infected", "MPXV clade IIa infected")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$Treatment <- relevel(dds$Treatment , ref = "MPXV clade I infected")
dds$Treatment


#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)


#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"normalized_counts_keratinocyte_IIa_vs_I.csv")
#####################################
# Boxplot

# Log2 transformation for count data
count_matrix <- counts(dds) + 1  # Adding 1 to avoid log(0)
log2_count_matrix <- log2(count_matrix)
boxplot(log2_count_matrix, outline = FALSE, main = "Boxplot of Log2-transformed Count Data")

# Log2 transformation for normalized count data
normalized_counts <- counts(dds, normalized = TRUE)
log2_normalized_counts <- log2(normalized_counts + 1)  # Adding 1 to avoid log(0)
boxplot(log2_normalized_counts, outline = FALSE, main = "Boxplot of Log2-transformed Normalized Count Data")
##

# Boxplot for count data
par(mfrow = c(1, 2))  # Create a 1x2 grid for side-by-side plots
boxplot(counts(dds), main = "Count Data", xlab = "Treatment", ylab = "Counts", col = "lightblue")

# Boxplot for normalized count data
boxplot(normalize_counts, main = "Normalized Count Data", xlab = "Treatment", ylab = "Normalized Counts", col = "lightgreen")

# Reset the layout
par(mfrow = c(1, 1))

####################################
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#for MPXV clade IIa infected
deseq_results_IIa <- results(dds, contrast=c("Treatment", "MPXV clade IIa infected", "MPXV clade I infected"))
deseq_results_IIa <- as.data.frame(deseq_results_IIa)
head(deseq_results_IIa)
deseq_results_IIa_0.05 <- results(dds,contrast=c("Treatment", "MPXV clade IIa infected", "MPXV clade I infected"), alpha = 0.05)
deseq_results_IIa_0.05 <- as.data.frame(deseq_results_IIa_0.05)
head(deseq_results_IIa_0.05)


# create histogram plot of p-values
hist(deseq_results_IIa_0.05$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "Frequencies of padj-values")


deseq_results <- deseq_results_IIa_0.05
summary(deseq_results)


#baseMean         log2FoldChange         lfcSE              stat               pvalue            padj      
#Min.   :    0.309   Min.   :-6.02053   Min.   :0.04332   Min.   :-22.29339   Min.   :0.0000   Min.   :0.000  
#1st Qu.:    2.056   1st Qu.:-0.42893   1st Qu.:0.36149   1st Qu.: -0.57022   1st Qu.:0.2559   1st Qu.:0.292  
#Median :   10.394   Median : 0.03209   Median :0.69160   Median :  0.02914   Median :0.5316   Median :0.673  
#Mean   :   57.049   Mean   : 0.06721   Mean   :1.22547   Mean   :  0.05383   Mean   :0.5144   Mean   :0.584  
#3rd Qu.:   35.265   3rd Qu.: 0.46585   3rd Qu.:1.69303   3rd Qu.:  0.68085   3rd Qu.:0.7862   3rd Qu.:0.884  
#ax.   :24141.945   Max.   :20.81238   Max.   :4.40736   Max.   : 25.96353   Max.   :1.0000   Max.   :1.000  
#NA's   :1        NA's   :10203


#Researchers are often interested in minimizing the number of false discoveries. 
#Change DEseq object to R object(dataframe)
deseq_results <- as.data.frame(deseq_results)
class(deseq_results)
head(deseq_results)


#Extract the most differentially expressed genes due to the treatment.
#Select genes with a significant change in gene expression (adjusted p-value below 0.05)
#And log2fold change <1 and >1

#Step 1: filter based on p adjusted value
filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)
head(filtered)
dim(filtered)
#Step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange)>1)
dim(deseq_results)
dim(filtered)

#make queries

#Save the deseg result. We will save the both the original data(res) and the filtered one(hits)
write.csv(deseq_results,"deseq_result_keratinocyte_IIa_vs_I.csv")
write.csv(filtered,"filtered_keratinocyte_IIa_vs_I.csv")



###################################################
#Batch effect
# Load necessary libraries
library(DESeq2)
BiocManager::install("sva")
BiocManager::install("svaCombat")
library(sva)
install.packages("svaComBat")
library(svaCombat)

# Define your metadata

sample_info2$batch <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
metadata <- sample_info2
# Load count matrix (make sure it's properly formatted)
# count_matrix <- read.table("your_count_matrix.txt", header = TRUE, row.names = 1)
batch_info <- sample_info2$batch  # Replace 'BatchColumn' with your actual batch column name
# Create a model matrix
model_matrix <- model.matrix(~ batch_info, data = sample_info2)

# Perform surrogate variable analysis (SVA)
sva_results <- sva(sample_info2, mod = model_matrix)
sva_results <- sva(sample_info2, mod = ~batch_info)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info2, design = ~ Treatment)

# Check for batch effects using principal component analysis (PCA)
vsd <- vst(dds)
vsd <- vst (dds,blind = FALSE)
pca <- plotPCA(vsd, intgroup = "batch")

vsd_matrix <- assay(vsd)
vsd_matrix
#vsd <- as.data.frame(vsd)
vsd
# Correct for batch effects using ComBat
combat_data <- ComBat(vsd_matrix, batch = sample_info2$batch)
boxplot(vsd_matrix)
boxplot(combat_data)
combat_data_rounded <- round(combat_data)
boxplot(combat_data_rounded) #outline = FALSE)
boxplot(count_data)
ComBat_seq(count_data,batch_info, group = Treatment, covar_mod = Cell_type,full_mod = ~Treatment+Cell_type,shrink = FALSE,
           shrink.disp = FALSE,gene.subset.n = NULL)
####################################
# Log2 transformation for count data
count_matrix <- counts(dds) + 1  # Adding 1 to avoid log(0)
log2_count_matrix <- log2(count_matrix)
boxplot(count_data, outline = FALSE, main = "Boxplot of Log2-transformed Count Data")

#cpm 
#deseq2
#combat
#######################
#set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
# Create a DESeqDataSet object using rounded combat_data
#set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Cell_type <- factor(sample_info$Cell_type)

# Convert non-integer values to integers in count data
#count_data <- round(count_data)
#head(count_data)

# Create a new count data object
combat_data_rounded <- as.matrix(combat_data_rounded)
head(combat_data_rounded)


unique(sample_info$Treatment)
unique(sample_info$Cell_type)





# Perform differential expression analysis with DESeq2 on the batch-corrected data
dds_combat <- DESeqDataSetFromMatrix(countData = combat_data_rounded, colData = sample_info2, design = ~ Treatment)
dds_combat <- DESeq(dds_combat)


#variance transformed values to generate a pca plot 
vsd <- vst (dds,blind = FALSE)

vsdrld <- rlog(dds,blind = FALSE)
head(assay(vsd),3)



#use transformed values to generate a pcs plot 
plotPCA(vsd,intgroup=c("Cell_type", "Treatment"))
# You can then extract differentially expressed genes, if needed.
# For example:
# res <- results(dds_combat)
# top_genes <- subset(res, padj < 0.05)

# Further downstream analysis can be performed with the results.

#########################################################

#visualization




#dispersion plot

plotDispEsts(dds)

#PCA
#pca stands for principal component analysis.
#it is a dimensionality reduction technique and in gene expression analysis, 139 #it can be used to explain the variance in gene expression datasets.
#to generate the pca plot we will first perform a variance stabilizing transformation.
#we will use the vst function in deseq.
#after that we use the transformed values to plot the PCA

#variance transformed values to generate a pca plot 
vsd <- vst (dds,blind = FALSE)

vsdrld <- rlog(dds,blind = FALSE)
head(assay(vsd),3)



#use transformed values to generate a pcs plot 
plotPCA(vsd,intgroup=c("Cell_type", "Treatment"))

#Heatmaps

#R Package: pheatmap
#Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts.

#generate the distance matrix
sampledists <- dist (t (assay(vsd)))
sampleDistMatrix <- as.matrix(sampledists)

# Define custom row and column names (replace these with your desired names)
#sample_info1 <- read.csv("meta_keratinocyte.csv")
custom_row_names <- paste(sample_info1$Sample.Name, sample_info1$Treatment, sep = "_")
custom_col_names <- paste(sample_info1$Sample.Name, sample_info1$Treatment, sep = "_")
rownames(sampleDistMatrix) <- custom_row_names
colnames(sampleDistMatrix) <- custom_col_names

#set color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Purples")))(255)#Blues,Purples
#generate the heatmap
pheatmap(sampleDistMatrix,clustering_distance_rows=sampledists,
         clustering_distance_cols = sampledists, col=colors)



#Heatmap of log transformed normalized counts. We will use the top 20 genes.
#Heatmap of Z scores. We will use the top 20 genes
# Heatmap of log transformed normalized counts . we will use the top 20 genes

#top 20 genes
top_hits <- deseq_results[order(deseq_results$padj),][1:20,]
top_hits <-row.names(top_hits)
top_hits
write.csv(top_hits,"top_hits_IIa_vs_I_20.csv")
class(top_hits)
# Order the values in the top_hits character vector
top_hits <- top_hits[order(top_hits)]

# Create the heatmap using the ordered top_hits
pheatmap(assay(rld)[top_hits, ], 
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE)

top_geneid <- read.csv("ensemble_geneid.csv")


# Filter top_geneid to only include rows where gene_id matches the top_hits
top_gene_names <- top_geneid[top_geneid$gene_id %in% top_hits, ]
top_gene_names$gene_name


# Save the filtered data to a CSV file
write.csv(top_gene_names, file = "top_gene_names_IIa_vs_I_20.csv", row.names = FALSE)


rld <- rlog(dds,blind = FALSE)

#Load the gene names from the CSV file
#top_gene_names <- read.csv("top_gene_names_IIa_vs_I_20.csv")

# Create a vector of gene names in the same order as the top_hits
gene_names_ordered <- top_gene_names$gene_name
show(gene_names_ordered)

# Create a heatmap with row names replaced by gene names
heatmap_data <- assay(rld)[top_hits, ]
ordered_gene_ids <- rownames(heatmap_data)[order(rownames(heatmap_data))]
heatmap_data_ordered <- heatmap_data[ordered_gene_ids, ]
colnames(heatmap_data_ordered) <- custom_col_names
pheatmap(heatmap_data_ordered, 
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE)

# Create the heatmap with ordered rows
rownames(heatmap_data_ordered) <- gene_names_ordered
pheatmap(heatmap_data_ordered, 
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE)

pheatmap(heatmap_data_ordered)

pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE, show_rownames = TRUE,cluster_cols = FALSE)
pheatmap(assay(rld)[top_hits,])

##############
annot_info <- as.data.frame(colData(dds)[, c("Cell_type", "Treatment")])
# Define colors for annotations (change these to your desired colors)
annotation_colors <- list(Cell_type = c("karatinocyte" = "red"),
                          Treatment = c("mock" = "green", "MPXV clade IIa infected" = "purple"))

# Use annotation_colors in pheatmap
pheatmap(heatmap_data_ordered, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE,
         annotation_col = annot_info, annotation_colors = annotation_colors)
################
annot_info <- as.data.frame(colData(dds)[,c("Cell_type", "Treatment")])
pheatmap(heatmap_data_ordered,cluster_rows = FALSE,show_rownames = TRUE,cluster_cols = FALSE,
         annotation_col = annot_info)


# Calculate Z-scores for the top 20 genes
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits, ]

# Create a heatmap of Z-scores with gene names as row names
rownames(zscore_subset) <- gene_names_ordered

pheatmap(zscore_subset)

#heatmap of Z scores . we will use the top 20 genes
cal_z_score <- function(x){(x-mean(x)) / sd(x)}


zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
zscore_subset <- zscore_all [top_hits,]
pheatmap(zscore_subset)

#MA plot
plotMA(dds,ylim = c (-2,2))

# Identify available coefficient names
coeff_names <- resultsNames(dds)

# Print the coefficient names
print(coeff_names)


#remove the noise 

#[4] "Treatment_MPXV.clade.IIa.infected_vs_MPXV.clade.I.infected" 

resLFC <- lfcShrink(dds, coef ="Treatment_MPXV.clade.IIa.infected_vs_MPXV.clade.I.infected" , type = "apeglm")

resplotMA(resLFC,ylim=c(-2,2))

#Volcano Plot

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)


# Create a Volcano Plot
ggplot(resLFC, aes(x = log2FoldChange, y = -log10(pvalue))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(pvalue) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(pvalue) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  
  # Customize plot labels and add the header
  labs(
    title = "HeatMap of Treatment_MPXV.clade.IIa.infected_vs_MPXV.clade.I.infected", # Add the header here
    x = "Log2 Fold Change",
    y = "-log10(P-Value)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

Upregulated <- resLFC[resLFC$log2FoldChange > 1 & resLFC$pvalue < 1.30103, ]
top_UP_hits1 <- Upregulated[order(Upregulated$pvalue),]
top_UP_hits <- Upregulated[order(Upregulated$pvalue),][1:20,]
top_UP_hits1 <-row.names(top_UP_hits1)
top_UP_hits <-row.names(top_UP_hits)
top_UP_hits
#write.csv(top_UP_hits,"top_UP_hits_IIa_vs_I_20.csv")
class(top_UP_hits)

top_UP_hits1 <- top_UP_hits1[order(top_UP_hits1)]
top_UP_hits <- top_UP_hits[order(top_UP_hits)]

top_geneid <- read.csv("ensemble_geneid.csv")


# Filter top_geneid to only include rows where gene_id matches the top_UP_hits
top_gene_names1 <- top_geneid[top_geneid$gene_id %in% top_UP_hits1, ]
top_gene_names <- top_geneid[top_geneid$gene_id %in% top_UP_hits, ]
top_gene_names$gene_name
top_gene_names1$gene_name

# Save the filtered data to a CSV file
write.csv(top_gene_names1, file = "UP.all_hits_gene_names_IIa_vs_I_20.csv", row.names = FALSE)
write.csv(top_gene_names, file = "UP.20_hits_gene_names_IIa_vs_I_20.csv", row.names = FALSE)


#up
#[1] "SERPINB1" "KLF6"     "ENO1"     "NDRG1"    "SLC2A1"   "PLAU"     "AHNAK"   
#[8] "BHLHE40"  "LDHA"     "TXN"      "IER3"     "FTH1"     "DDIT4"    "SFN"     
#[15] "PLEC"     "AKR1C1"   "S100A16"  "LAMB3"    "TXNRD1"   "SRXN1" 

#top hits
#" "EGR1"      "CLU"       "HIST1H2BJ" "KRT17"     "DNAJA4"    "MCL1"      "DDIT4"     "CXCL8"    
#[9] "SPRR1B"    "FOS"       "JUN"       "KRT6B"     "KRT5"      "KRT14"     "HIST1H4C"  "NRARP"    
#[17] "KRT6A"     "HLA-A"     "HLA-B"     "HIST1H2BG"

#down
# [1] "EGR1"      "CLU"       "HIST1H2BJ" "KRT17"     "DNAJA4"    "MCL1"      "HNRNPDL"  
#[8] "CXCL8"     "SPRR1B"    "FOS"       "JUN"       "KRT6B"     "KRT5"      "KRT14"    
#[15] "HIST1H4C"  "NRARP"     "KRT6A"     "HLA-A"     "HLA-B"     "HIST1H2BG"
Downregulated <- resLFC[resLFC$log2FoldChange < -1 & resLFC$pvalue < 1.30103, ]
top_hits1 <- Downregulated[order(Downregulated$pvalue),]
top_hits <- Downregulated[order(Downregulated$pvalue),][1:20,]
top_hits1 <-row.names(top_hits1)
top_hits <-row.names(top_hits)
top_hits
#write.csv(top_hits,"test.top_hits_IIa_vs_I_20.csv")
class(top_hits)
top_hits1 <- top_hits1[order(top_hits1)]
top_hits <- top_hits[order(top_hits)]

top_geneid <- read.csv("ensemble_geneid.csv")


# Filter top_geneid to only include rows where gene_id matches the top_hits
top_gene_names1 <- top_geneid[top_geneid$gene_id %in% top_hits1, ]
top_gene_names <- top_geneid[top_geneid$gene_id %in% top_hits, ]
top_gene_names1$gene_name
top_gene_names$gene_name


# Save the filtered data to a CSV file
write.csv(top_gene_names1, file = "DOWN.all_gene_names_IIa_vs_I_20.csv", row.names = FALSE)
write.csv(top_gene_names, file = "DOWN.20_gene_names_IIa_vs_I_20.csv", row.names = FALSE)

