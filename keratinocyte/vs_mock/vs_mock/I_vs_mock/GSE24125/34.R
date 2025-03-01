library(GEOquery)
library(limma)
library(umap)

gset <- getGEO("GSE11234", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6764", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("XXXXXXXXXXXXXXX1111111122222XXXX1000001XXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXX")

sml <- strsplit(gsms, split="")[[1]]

#boxplot(ex)

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

########################################################################
normalized.data <- rma(ex)
########################################################################
# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Mpox","Mock","Killed"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c("Mpox-Mock", "Mpox-Killed", "Killed-Mock")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# Order samples by group
ord <- order(gs)

# Define a palette for the boxplot colors
my_palette <- c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
                "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666")

# Set up the plotting environment
par(mar=c(7, 4, 2, 1))

# Create the boxplot
title <- paste("Boxplot based on infection type ", "/", annotation(gset), sep="")
boxplot(ex[, ord], boxwex=0.6, notch=TRUE, main=title, outline=FALSE, las=2, col=my_palette[gs[ord]])

####################################################################################################
ex1 <- normalizeQuantiles(ex)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

########################################################################
########################################################################
# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Mpox","Mock","Killed"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c("Mpox-Mock", "Mpox-Killed", "Killed-Mock")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# Order samples by group
ord <- order(gs)

# Define a palette for the boxplot colors
my_palette <- c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
                "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666")

# Set up the plotting environment
par(mar=c(7, 4, 2, 1))

# Create the boxplot
title <- paste("Normalized Boxplot based on infection type ", "/", annotation(gset), sep="")
boxplot(ex[, ord], boxwex=0.6, notch=TRUE, main=title, outline=FALSE, las=2, col=my_palette[gs[ord]])

# assumption is that most genes are not differentially expressed.
hist(tT$P.Value , col = "#9c0b6e", border = "white", xlab = "P-Value",
     ylab = "Number of genes", main = "P Value distribution")

hist(tT$adj.P.Val , col = "#9c0b6e", border = "white", xlab = "Adjusted P value",
     ylab = "Number of genes", main = "P adjusted Value distribution")
