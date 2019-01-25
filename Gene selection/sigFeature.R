#######This algorithm could take huge time
library(sigFeature)
####################################

x <- as.matrix(TCGALUAD[,3:14479])
y <- as.matrix(TCGALUAD$status)

sigrank <- sigFeature(x,y)
