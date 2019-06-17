###For GSEA, we first calculated each gene's log foldchange using edgeR package by comparing high-risk to low-risk group.
###Then, using clusterprofiler to perform GSEA with Reactome database.
###The Reactome gmt file can be downloaded from http://software.broadinstitute.org/gsea/msigdb/

library(clusterProfiler)
library(enrichplot)

df <- gseainput
df
df.id<-bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(df.id)
easy.df<-merge(df,df.id,by="SYMBOL",all=F)

sortdf<-easy.df[order(easy.df$foldChange, decreasing = T),]
head(sortdf)
gene.expr = sortdf$foldChange
names(gene.expr) <- sortdf$ENTREZID
head(gene.expr)

genesets <- read.gmt("c2.cp.reactome.v6.2.symbols.gmt.txt")

kk <- GSEA(geneList,TERM2GENE = genesets,
           exponent = 1, 
           nPerm = 1000, 
           minGSSize = 10,
           maxGSSize = 500, 
           pvalueCutoff = 1, 
           verbose = TRUE, 
           seed = TRUE,
           by = "fgsea")

kk
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]

gseaplot2(kk, row.names(sortkk)[1:4])
ggsave('1.pdf')

gseaplot2(kk, row.names(sortkk)[584:587])
ggsave('2.pdf')
gseaplot2

sortkk <- filter(sortkk,sortkk$pvalue<0.05)
write.csv(sortkk,'gsea.csv')




