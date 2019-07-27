library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(caret)

#######download data
gset <- getGEO("GSE72094", GSEMatrix =TRUE, getGPL = TRUE, AnnotGPL = TRUE
)
if (length(gset) > 1) idx <- grep("GPL15048", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
exprdf<-data.frame(exprs(gset))
head(exprdf)

exprdf [1:4,1:4]


####prepare clinical data
clinicaldata <- pData(gset)

clinicaldata <- clinicaldata[,c(2,14,15,19,20,21)]
clinicaldata$times <- str_split(clinicaldata$characteristics_ch1.10,':',simplify = T)[,2]
clinicaldata$times <- as.numeric(clinicaldata$times)
clinicaldata_new <- na.omit(clinicaldata)
clinicaldata_new$characteristics_ch1.9 <- as.character(clinicaldata_new$characteristics_ch1.9)
clinicaldata_new$status <- ifelse(clinicaldata_new$characteristics_ch1.9 %like% 'Alive',0,1)

############extract genes
##check probe.xlsx
probe$probe
aaa <- intersect(probe$probe,colnames(exprdf))

expr_need <- exprdf[,aaa]

expr_need <- as.data.frame(t(expr_need))
expr_need <- expr_need[,clinicaldata_new$geo_accession]
expr_need$probe <- row.names(expr_need)
expr_need_gene <- merge(probe,expr_need ,by='probe')
expr_need_gene <- expr_need_gene [,-1]
expr_need_gene_expression <- aggregate(.~Gene,expr_need_gene,mean)####average gene expression for multiple probes
row.names(expr_need_gene_expression) <- expr_need_gene_expression$Gene
expr_need_gene_expression <- expr_need_gene_expression[,-1]


#####calculate riskscore


lasso_index <- fread('lasso.txt') 

lasso_index <- as.data.frame(lasso_index)
row.names(lasso_index) <- lasso_index[,1]

expr_need_gene_expression <- expr_need_gene_expression[lasso_index$Gene,]

lasso_index <- lasso_index[,-1]


expr_need_gene_expression <- as.data.frame(t(expr_need_gene_expression))
expr_need_gene_expression <- as.matrix(expr_need_gene_expression)
riskscore <- expr_need_gene_expression %*% lasso_index 


riskscore <- as.data.frame(riskscore)
colnames(riskscore) <- 'riskscore'
riskscore$sampleid <- row.names(riskscore)


colnames(clinicaldata_new)
riskscore$geo_accession <- row.names(riskscore)



expr_need_gene_expression_risk <- merge(clinicaldata_new,riskscore,by='geo_accession')
expr_need_gene_expression_risk <- filter(expr_need_gene_expression_risk,expr_need_gene_expression_risk$times>=30)
#######calculate c-index
rcorr.cens(-expr_need_gene_expression_risk$riskscore,Surv(expr_need_gene_expression_risk$times,expr_need_gene_expression_risk$status))[[1]]


#######grouping
res.cutsur <- surv_cutpoint(expr_need_gene_expression_risk , time = "times", 
                            event = "status", 
                            variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]

expr_need_gene_expression_risk$risk<-as.vector(ifelse(expr_need_gene_expression_risk$riskscore>=res.cut,"high","low"))

sum(expr_need_gene_expression_risk$risk=="high")
sum(expr_need_gene_expression_risk$risk=="low")

expr_need_gene_expression_risk$risk <- as.factor(expr_need_gene_expression_risk$risk)

expr_need_gene_expression_risk$risk <- relevel(expr_need_gene_expression_risk$risk,ref = 'low')
####cox analysis
cox <- coxph(Surv(times, status) ~ risk, data =expr_need_gene_expression_risk)
coxSummary = summary(cox)
#check results
coxSummary



sum(expr_need_gene_expression_risk$risk=="high")
sum(expr_need_gene_expression_risk$risk=="low")

#draw survival curve
Sur <- Surv(expr_need_gene_expression_risk$time,expr_need_gene_expression_risk$status)
sfit <- survfit(Sur ~ risk,data=expr_need_gene_expression_risk)
group <- factor(expr_need_gene_expression_risk$risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit,
           conf.int=F, 
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=155)","Low risk (n=231)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'GSE72094 (n=386)')
ggsave('GSE72094.pdf',width = 7,height = 7)

###multivariate cox

expr_need_gene_expression_risk$age <- as.numeric(str_split(expr_need_gene_expression_risk$characteristics_ch1.5,':',simplify = T)[,2])


expr_need_gene_expression_risk$gender <- as.factor(as.character(str_split(expr_need_gene_expression_risk$characteristics_ch1.4,':',simplify = T)[,2]))


expr_need_gene_expression_risk$stage <- as.factor(substr(expr_need_gene_expression_risk$characteristics_ch1.11,7,8))




mucox <- coxph(Surv(times, status) ~ risk+age+gender+stage, data =expr_need_gene_expression_risk)
coxSummarymu = summary(mucox)
#check results
coxSummarymu



