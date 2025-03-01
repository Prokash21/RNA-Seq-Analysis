x <- count_data[count_data$]
x <- count_data1[count_data1$gene_id %in% top_hits, ]

write.csv(x,file = "complexheat")
heatmap(x, Rowv = NULL, Colv = if(symm)"Rowv" else NULL,
        distfun = dist, hclustfun = hclust,
        reorderfun = function(d, w) reorder(d, w),
        add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
        scale = c("row", "column", "none"), na.rm = TRUE,
        margins = c(5, 5), ColSideColors, RowSideColors,
        cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
        labRow = NULL, labCol = NULL, main = NULL,
        xlab = NULL, ylab = NULL,
        keep.dendro = FALSE, verbose = getOption("verbose"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
install.packages("ComplexHeatmap")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

x1<-read.csv("complexheat",header=TRUE)
x <- x[, -1]
rownames(x) <- x[, 1]
x <- x[, -1]
x <- as.matrix(x)

Heatmap3D(x)
Heatmap(x)
HeatmapAnnotation(x)
HeatmapList(x)
heat