if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
res.df <- as.data.frame(resLFC)
res.df$symbol <- mapIds(org.Hs.eg.db,keys = rownames(res.df),keytype= "ENSEMBL",column=
                          "SYMBOL")
res.df

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

install.packages("textshaping")
library(EnhancedVolcano)


EnhancedVolcano(res.df, x = "log2FoldChange" , y = "padj" ,lab = res.df$symbol)
write.csv(res.df, file = "Graphprism.csv" )
write.csv()

EnhancedVolcano(res.df, x = "log2FoldChange" , y = "padj" ,lab = res.df$symbol,
                pCutoff= 0.05 , FCcutoff = 1)


select= c("DNAJA1" , "PPP1R15A" , "HSP90AB1" , "DDX5", "HSPA8"  ,"PKP1"  ,   "GSTP1" , 
          "SCD"  ,    "GAPDH"  ,  "KRT17" ,   "FGFBP1"   ,"SH3BGRL3" )

EnhancedVolcano(res.df ,x = "log2FoldChange" , y = "padj" ,lab = res.df$symbol,
                pCutoff= 0.05 , FCcutoff = 1,selectLab = select) 



write.csv(res.df,"DE.res.df.result.csv")



