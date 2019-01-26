#######This algorithm could take huge time(hours to days)
library(sigFeature)
####################################
x <- as.matrix(TCGALUAD[,3:ncol(TCGALUAD)])
y <- as.matrix(TCGALUAD$status)
sigrank <- sigFeature(x,y)
###################### top1000 variable extraction
signame_top1000 <- colnames(x[,sigrank[1:1000]])
