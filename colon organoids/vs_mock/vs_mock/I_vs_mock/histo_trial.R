#############################################################

count_data1<- read.csv("normalized_counts_keratinocyte_I_vs_mock.csv", header = TRUE)

#top_genes_series<- top_gene_names$gene_id
#x <- count_data1[count_data1$gene_id %in% top_hits, ]
top_genes_count<-count_data1 [count_data1$gene_id %in% top_hits,]



write.csv(top_genes_count,"count_of_top_genes.csv")

## Create the heatmap using the ordered top_hits
#heatmap
x<-as.data.frame(assay(rld)[top_hits, ]) 

write.csv(x,file="count_of_rld.csv")

#histogram
y<-as.data.frame(deseq_results_IIa_0.05$padj) 
write.csv(y,file="histogram_padj_data.csv")

ss<-hist(deseq_results_IIa_0.05$padj)
hist(deseq_results_IIa_0.05$pad)
