library(tidyverse)
library(survival)
##################################
TCGACoxoutput=data.frame()
for(i in colnames(TCGALUAD[,3:ncol(TCGALUAD)])){
  print(paste("trying for",i,"gene"))
  cox <- coxph(Surv(times, status) ~ TCGALUAD[,i], data = TCGALUAD)
  coxSummary = summary(cox)
  TCGACoxoutput=rbind(TCGACoxoutput,cbind(gene=i,HR=round(coxSummary$coefficients[,"exp(coef)"],3),
                                          z=coxSummary$coefficients[,"z"],
                                          pvalue=round(coxSummary$coefficients[,"Pr(>|z|)"],5),
                                          lower95CI=round(coxSummary$conf.int[,3],3),
                                          upper95CI=round(coxSummary$conf.int[,4],3)))
}
for(i in c(2:6)){
  TCGACoxoutput[,i] <- as.numeric(TCGACoxoutput[,i])
}
TCGACoxoutput0.005 <- arrange(TCGACoxoutput,pvalue)  %>% 
  filter(pvalue < 0.005) 
