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
dds$Treatment <- relevel(dds$Treatment , ref = "MPXV clade I infected")
dds$Treatment


#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#for MPXV clade IIa infected
deseq_results_IIa <- results(dds, contrast=c("Treatment", "MPXV clade IIa infected", "MPXV clade I infected"))
deseq_results_IIa <- as.data.frame(deseq_results_IIa)
head(deseq_results_IIa)
deseq_results_IIa_0.05 <- results(dds,contrast=c("Treatment", "MPXV clade IIa infected", "MPXV clade I infected"), alpha = 0.05)
deseq_results_IIa_0.05 <- as.data.frame(deseq_results_IIa_0.05)
head(deseq_results_IIa_0.05)



#dds$group <- factor(paste0(dds$Treatment))
#design(dds) <- ~ group
#dds <- DESeq(dds)
#resultsNames(dds)
#results(dds,contrast=c("Treatment", "MPXV clade IIa infected", "MPXV clade I infected") )
res <- results(dds)
metadata(res)$alpha
metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


#Count outlier detection
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

# create histogram plot of p-values
hist(deseq_results_IIa_0.05$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "Frequencies of padj-values")


deseq_results <- deseq_results_IIa_0.05

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

#Step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange)>1)
dim(deseq_results)
dim(filtered)

#make queries

#Save the deseg result. We will save the both the original data(res) and the filtered one(hits)
write.csv(deseq_results,"deseq_result_keratinocyte_IIa_vs_I.csv")
write.csv(filtered,"filtered_keratinocyte_IIa_vs_I.csv")

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"normalized_counts_keratinocyte_IIa_vs_I.csv")

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

rld <- rlog(dds,blind = FALSE)
head(assay(vsd),3)

##########################################################################
# this gives log2(n + 1)
#BiocManager::install("vsn")
ntd <- normTransform(dds)
library("vsn")

#meanSdPlot(assay(ntd))
#meanSdPlot(assay(vsd))
#meanSdPlot(assay(rld))


pcaData <- plotPCA(vsd, intgroup=c("Cell_type", "Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Cell_type, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#####################################################################################

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

# Load the gene names from the CSV file
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

####################################################################
pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE, show_rownames = TRUE,cluster_cols = FALSE)
pheatmap(assay(rld)[top_hits,])


annot_info <- as.data.frame(colData(dds)[,c("Cell_type", "Treatment")])
pheatmap(heatmap_data_ordered,cluster_rows = FALSE,show_rownames = TRUE,cluster_cols = FALSE,
         annotation_col = annot_info)

##############
# Calculate Z-scores for the top 20 genes
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits, ]

# Create a heatmap of Z-scores with gene names as row names
rownames(zscore_subset) <- gene_names_ordered

pheatmap(zscore_subset)
##################
#heatmap of Z scores . we will use the top 20 genes
cal_z_score <- function(x){(x-mean(x)) / sd(x)}


zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
zscore_subset <- zscore_all [top_hits,]
pheatmap(zscore_subset)




#MA plot
plotMA(dds,ylim = c (-2,2))

#remove the noise 
#resLFC <- lfcShrink(dds ,coef = "Treatment_mock_vs_MPXV clade I infected", type="apeglm")

##################
# Identify available coefficient names
coeff_names <- resultsNames(dds)

# Print the coefficient names
print(coeff_names)

#[1] "Intercept"                                                                        
#[2] "Cell_line_Normal.Human.Epidermal.Keratinocytes_vs_iPS.cell.derived.colon.organoids"
#[3] "Treatment_MPXV.clade.I.infected_vs_mock"                                           
#[4] "Treatment_MPXV.clade.IIb.infected_vs_mock"                                         
#[5] "Treatment_MPXV.clade.IIa.infected_vs_mock"                                         
#    "Treatment_MPXV.clade.I.infected_vs_MPXV.clade.IIb.infected"

#[3] "Treatment_MPXV.clade.I.infected_vs_mock"   

resLFC <- lfcShrink(dds, coef ="Treatment_MPXV.clade.IIa.infected_vs_MPXV.clade.I.infected" , type = "apeglm")

plotMA(resLFC,ylim=c(-2,2))

########################


#Volcano Plot

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
#label the genes 


resLFC$diffexpressed <- "NO"

resLFC$diffexpressed[resLFC$log2FoldChange> 1.0 & resLFC$padj<0.05]<- "UP"

resLFC$diffexpressed[resLFC$log2FoldChange< -1.0 & resLFC$padj<0.05]<- "DOWN"

resLFC$delabel<-NA


#ggplot(data=resLFC,aes(x=log2FoldChange,y =-log10(pvalue),col=diffexpressed,label=delabel))+
#  geom_point()+
#  theme_minimal()+
#  geom_text_repel()+
#  scale_color_manual(values=c("blue","black","red"))+
#  theme(text = element_text(size = 20))
#########################################



#####################################################
# Create a Volcano Plot
ggplot(resLFC, aes(x = log2FoldChange, y = -log10(pvalue))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange >= 1.0 & -log10(pvalue) >= 1.3, "Upregulated",
                                ifelse(log2FoldChange <= -1.0 & -log10(pvalue) >= 1.3, "Downregulated", "Not Significant"))),
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

#ggplot(resLFC,aes(x= log2FoldChange , y = -log10(padj), color=diffexpressed)) +
#  geom_point()
#####################################################
#ggsave(file.path(output_path, "volcano_plot.png"), plot = volcano_plot)

# Order the dataframe for the top 20 upregulated genes
ordered_up_resLFC <- resLFC[resLFC$diffexpressed == "UP" & resLFC$log2FoldChange > 1 & resLFC$padj < 0.05, ]
ordered_up_resLFC <- ordered_up_resLFC[order(ordered_up_resLFC$diffexpressed), ][1:20, ]

ordered_up_resLFC <-row.names(ordered_up_resLFC)
ordered_up_resLFC
write.csv(ordered_up_resLFC, file = "IIa_vs_I.ordered_up_resLFC.csv", row.names = FALSE)

top_geneid <- read.csv("ensemble_geneid.csv")


# Filter top_geneid to only include rows where gene_id matches the top_hits
top_gene_names <- top_geneid[top_geneid$gene_id %in% ordered_up_resLFC, ]
top_gene_names$gene_name


# Save the filtered data to a CSV file
write.csv(top_gene_names, file = "test.UP_gene_names_IIa_vs_I_20.csv", row.names = FALSE)

#write.csv(ordered_up_resLFC, file = "I.log2_ordered_up_resLFC.csv", row.names = FALSE)

# Order the dataframe for the top 20 downregulated genes
ordered_down_resLFC <- resLFC[resLFC$diffexpressed == "DOWN" & resLFC$log2FoldChange < -1 & resLFC$padj < 0.05, ]
ordered_down_resLFC <- ordered_down_resLFC[order(ordered_down_resLFC$diffexpressed), ][1:20, ]
ordered_down_resLFC <- row.names(ordered_down_resLFC)
ordered_down_resLFC
write.csv(ordered_down_resLFC, file = "IIa_vs_I.ordered_down_resLFC.csv", row.names = FALSE)
#write.csv(ordered_down_resLFC, file = "I.log2_ordered_down_resLFC.csv", row.names = FALSE)

# Filter top_geneid to only include rows where gene_id matches the top_hits
top_gene_names <- top_geneid[top_geneid$gene_id %in% ordered_down_resLFC, ]
top_gene_names$gene_name


# Save the filtered data to a CSV file
write.csv(top_gene_names, file = "DOWN_gene_names_IIa_vs_I_20.csv", row.names = FALSE)

