

# Install and load the VennDiagram package (if not already installed)
install.packages("VennDiagram")
library(VennDiagram)

# Define the gene sets
up_genes_I <- c(
  "ENSG00000001630", "ENSG00000008394", "ENSG00000034510", "ENSG00000044574", "ENSG00000067064",
  "ENSG00000067082", "ENSG00000080824", "ENSG00000081051", "ENSG00000086061", "ENSG00000091513",
  "ENSG00000099194", "ENSG00000107742", "ENSG00000108242", "ENSG00000109072", "ENSG00000109321",
  "ENSG00000109472", "ENSG00000109814", "ENSG00000109971", "ENSG00000110244", "ENSG00000110245"
)

up_genes_IIa <- c(
  "ENSG00000109321", "ENSG00000118523", "ENSG00000120129", "ENSG00000120738", "ENSG00000128016",
  "ENSG00000132002", "ENSG00000170345", "ENSG00000173110", "ENSG00000204389", "ENSG00000282885"
)

# Create a Venn diagram object
venn_obj <- venn.diagram(
  x = list(CladeI = up_genes_I, CladeIIa = up_genes_IIa),
  category.names = c("Clade I", "Clade IIa"),
  filename = NULL  # Use NULL to display the plot in the R console
)

# Display the Venn diagram
grid.draw(venn_obj)


# Define gene symbols for each treatment comparison
up_genes_I <- c(
  "ENSG00000001630", "ENSG00000008394", "ENSG00000034510", "ENSG00000044574", "ENSG00000067064",
  "ENSG00000067082", "ENSG00000080824", "ENSG00000081051", "ENSG00000086061", "ENSG00000091513",
  "ENSG00000099194", "ENSG00000107742", "ENSG00000108242", "ENSG00000109072", "ENSG00000109321",
  "ENSG00000109472", "ENSG00000109814", "ENSG00000109971", "ENSG00000110244", "ENSG00000110245"
)

# Define gene symbols for Clade IIb (if available)
# Leave it empty if there are no genes for Clade IIb

# Define gene symbols for Clade IIa
up_genes_IIa <- c(
  "ENSG00000109321", "ENSG00000118523", "ENSG00000120129", "ENSG00000120738", "ENSG00000128016",
  "ENSG00000132002", "ENSG00000170345", "ENSG00000173110", "ENSG00000204389", "ENSG00000282885"
)

# Create a list of gene sets for each treatment comparison
gene_sets <- list(
  up_genes_I = up_genes_I,
  up_genes_IIb = character(0),  # Empty list for Clade IIb if no genes available
  up_genes_IIa = up_genes_IIa
)


# Create a list of gene sets for each treatment comparison
gene_sets <- list(
  up_genes_I = c(
    "ENSG00000001630", "ENSG00000008394", "ENSG00000034510", "ENSG00000044574", "ENSG00000067064",
    "ENSG00000067082", "ENSG00000080824", "ENSG00000081051", "ENSG00000086061", "ENSG00000091513",
    "ENSG00000099194", "ENSG00000107742", "ENSG00000108242", "ENSG00000109072", "ENSG00000109321",
    "ENSG00000109472", "ENSG00000109814", "ENSG00000109971", "ENSG00000110244", "ENSG00000110245"
  ),
  up_genes_IIb = c(),  # Add gene symbols for Clade IIb here
  up_genes_IIa = c("ENSG00000109321", "ENSG00000118523", "ENSG00000120129", "ENSG00000120738", "ENSG00000128016",
                   "ENSG00000132002", "ENSG00000170345", "ENSG00000173110", "ENSG00000204389", "ENSG00000282885")   # Add gene symbols for Clade IIa here
)

# Create a list of gene sets for each treatment comparison
gene_sets <- list(
  up_genes_I = ordered_up_resLFC$GeneSymbol[ordered_up_resLFC$diffexpressed == "UP" & 
                                              ordered_up_resLFC$log2FoldChange > 0.1 & ordered_up_resLFC$padj < 0.05],
  up_genes_IIb = ordered_up_resLFC$GeneSymbol[ordered_up_resLFC$diffexpressed == "UP" & 
                                                ordered_up_resLFC$log2FoldChange > 0.1 & ordered_up_resLFC$padj < 0.05],
  up_genes_IIa = ordered_up_resLFC$GeneSymbol[ordered_up_resLFC$diffexpressed == "UP" & 
                                                ordered_up_resLFC$log2FoldChange > 0.1 & ordered_up_resLFC$padj < 0.05]
)

install.packages("VennDiagram")
library(VennDiagram)
# Create a Venn diagram object
venn.diagram(
  x = gene_sets,
  category.names = c("Clade I", "Clade IIa"),
  filename = "venn_diagram.png"  # You can specify a filename for the plot
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Add another visualization plot (e.g., a barplot of gene counts for each treatment)
# Here's an example using ggplot2, but you can customize this plot as needed
library(ggplot2)

# Create a dataframe for gene counts in each treatment
gene_counts <- data.frame(
  Treatment = c("Clade I", "Clade IIb", "Clade IIa"),
  GeneCount = c(
    length(gene_sets$up_genes_I),
    length(gene_sets$up_genes_IIb),
    length(gene_sets$up_genes_IIa)
  )
)

# Create a barplot
bar_plot <- ggplot(gene_counts, aes(x = Treatment, y = GeneCount)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(
    title = "Number of Upregulated Genes in Each Treatment",
    x = "Treatment",
    y = "Number of Genes"
  ) +
  theme_minimal()

# Save the barplot as an image (e.g., PNG)
ggsave("gene_count_barplot.png", plot = bar_plot, width = 6, height = 4)


