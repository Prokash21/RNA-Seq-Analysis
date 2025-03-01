# Load required libraries
library(DESeq2)
library(VennDiagram)


#set the working directory
setwd("D:/DEG/GSE219036_6_june/keratinocyte/vs_mock/venn_kera")


# Load your DESeq2 result files
cladeI_vs_mock <- read.csv("UP.all_hits_gene_names_I_vs_mock.csv")
cladeIIa_vs_mock <- read.csv("UP.all_hits_gene_names_IIa_vs_mock.csv")
cladeIIb_vs_mock <- read.csv("UP.all_hits_gene_names_IIb_vs_mock_20.csv")



# Extract upregulated genes from each comparison
upregulated_genes_cladeI <- cladeI_vs_mock$gene_name
upregulated_genes_cladeIIa <- cladeIIa_vs_mock$gene_name
upregulated_genes_cladeIIb <- cladeIIb_vs_mock$gene_name


# Create a list of upregulated genes from each comparison
upregulated_gene_lists <- list(CladeI = upregulated_genes_cladeI,
                               CladeIIa = upregulated_genes_cladeIIa,
                               cladeIIb = upregulated_genes_cladeIIb)

# Define colors for the Venn diagram
colors <- c("dodgerblue", "darkorange", "forestgreen")


# Create a Venn diagram
venn.plot <- venn.diagram(
  x = upregulated_gene_lists,
  category.names = c("I_vs_mock", "IIa_vs_mock", "IIb_vs_mock"),
  filename = NULL,
  col = colors  
)
library(gplots)
venn(upregulated_gene_lists)

# Plot the Venn diagram
pdf("venn_diagram.pdf")
grid.draw(venn.plot)
dev.off()


# Find common genes among the three comparisons
common_genes <- Reduce(intersect, list(upregulated_genes_cladeI, upregulated_genes_cladeIIa, upregulated_genes_cladeIIb))
# Create a data frame with common genes
common_genes_df <- data.frame(Gene_Name = common_genes)

# Write the common genes to a CSV file
write.csv(common_genes_df, "common_genes_keratinocyte.csv", row.names = FALSE)



