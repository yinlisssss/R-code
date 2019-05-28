##Use gender as example

tcgaclinical <- read.csv("TCGAclinical.csv",check.names = F)
ttrain<-TCGALUAD_LASSO
ttrain$riskscore <- best.trrainrs
res.cutsur <- surv_cutpoint(ttrain, time = "times", 
                            event = "status", 
                            variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]
risk<-as.vector(ifelse(best.trrainrs>=res.cut,"high","low"))
ttrain$risk<-risk
strtcga_age <- ttrain[,c(1,2,575)]
tcgaclinical_age <- tcgaclinical[,c(1,12)]

strtcga_age$X_PATIENT <- row.names(strtcga_age)
str_age_TCGA <- merge(strtcga_age,tcgaclinical_age,by="X_PATIENT")

MALE_SUR <- filter(str_age_TCGA,str_age_TCGA$gender.x=='MALE')
sum(MALE_SUR$risk=='high')
sum(MALE_SUR$risk=='low')
Sur_male <- Surv(MALE_SUR$times/365,MALE_SUR$status)
sfit_male <- survfit(Sur_male ~ risk,data=MALE_SUR )
group <- factor(MALE_SUR$risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur_male ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit_male,
           conf.int=F, 
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=42)","Low risk (n=187)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD male patients (n=229)')
