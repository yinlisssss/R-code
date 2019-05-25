

###########分层分析
#TCGA STAGE1
tcgaclinical <- read.csv("TCGAclinical.csv",check.names = F)
strtcga_stage <- ttrain[,c(1,2,575)]
tcgaclinical_stage <- tcgaclinical[,c(1,15)]

strtcga_stage$X_PATIENT <- row.names(strtcga_stage)
str_STAGE_TCGA <- merge(strtcga_stage,tcgaclinical_stage,by="X_PATIENT")

smoking_SUR <- filter(str_STAGE_TCGA,str_STAGE_TCGA$stage=='Stage I')
sum(smoking_SUR$risk=='high')
sum(smoking_SUR$risk=='low')
Sur_stage1 <- Surv(smoking_SUR$times/365,smoking_SUR$status)
sfit_stage1 <- survfit(Sur_stage1 ~ risk,data=smoking_SUR )
group <- factor(smoking_SUR $risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur_stage1 ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit_stage1,
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=26)","Low risk (n=238)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD stage I patients (n=264)')
ggsave('stage1tcga.pdf',width = 7,height = 7)


smoking_SUR <- filter(str_STAGE_TCGA,str_STAGE_TCGA$stage=='Stage II')
sum(smoking_SUR$risk=='high')
sum(smoking_SUR$risk=='low')
Sur_stage1 <- Surv(smoking_SUR$times/365,smoking_SUR$status)
sfit_stage2 <- survfit(Sur_stage1 ~ risk,data=smoking_SUR )
group <- factor(smoking_SUR $risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur_stage1 ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit_stage2 ,
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=20)","Low risk (n=96)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD stage II patients (n=116)')
ggsave('stage2tcga.pdf',width = 7,height = 7)


smoking_SUR <- filter(str_STAGE_TCGA,str_STAGE_TCGA$stage=='Stage III')
sum(smoking_SUR$risk=='high')
sum(smoking_SUR$risk=='low')
Sur_stage1 <- Surv(smoking_SUR$times/365,smoking_SUR$status)
sfit_stage3 <- survfit(Sur_stage1 ~ risk,data=smoking_SUR )
group <- factor(smoking_SUR $risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur_stage1 ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit_stage3 ,
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=26)","Low risk (n=53)"),
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD stage III patients (n=79)')
ggsave('stage3tcga.pdf',width = 7,height = 7)


###########TCGA GENDER

tcgaclinical <- read.csv("TCGAclinical.csv",check.names = F)
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
           conf.int=F, #置信区间
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
ggsave('TCGAmale.pdf',height = 7,width = 7)



MALE_SUR <- filter(str_age_TCGA,str_age_TCGA$gender.x=='FEMALE')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=37)","Low risk (n=226)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD female patients (n=263)')
ggsave('TCGAfemale.pdf',height = 7,width = 7)

37+42

#########AGE

strtcga_age <- ttrain[,c(1,2,575)]
tcgaclinical_age <- tcgaclinical[,c(1,11)]

strtcga_age$X_PATIENT <- row.names(strtcga_age)
str_age_TCGA <- merge(strtcga_age,tcgaclinical_age,by="X_PATIENT")
str_age_TCGA$AGE <- ifelse(str_age_TCGA$age_at_initial_pathologic_diagnosis.x>65,'Yes','No')


MALE_SUR <- filter(str_age_TCGA,str_age_TCGA$AGE=='No')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=35)","Low risk (n=196)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD <=65 patients (n=231)')
ggsave('TCGAage<=65.pdf',height = 7,width = 7)


###################TCGA其他
###############测试所有
ttest1_str1 <- ttest1[,c(1,2,575)]
ttest2_str2 <- ttest2[,c(1,2,575)]
allltest_str <- rbind(ttest1_str1,ttest2_str2)
clinical19188 <- read.csv("clinical19188.csv")
clinical19188 <- clinical19188[,c("Geo_accession","Gender")]
clinical30219 <- read.csv("clinical30219.csv")
clinical30219 <- clinical30219[,c("Geo_accession","Gender")]
clinical31210 <- read.csv("clinical31210.csv")
clinical31210 <- clinical31210[,c("Geo_accession","Gender")]
clinical37735 <- read.csv("clinical37745.csv")
clinical37735 <- clinical37735[,c("Geo_accession","Gender")]
clinical50081 <- read.csv("clinical50081.csv")
clinical50081 <- clinical50081[,c("Geo_accession","Gender")]
clinicalall <- rbind(clinical19188,clinical30219,clinical31210,clinical37735,clinical50081)
allltest_str$Geo_accession <- row.names(allltest_str)
geoallgender <- merge(allltest_str,clinicalall,by='Geo_accession')

#######male
MALE_SUR <- filter(geoallgender,geoallgender$Gender=='Male')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=134)","Low risk (n=169)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'GEO male patients (n=303)')
ggsave('GEOmale.pdf',width = 7,height = 7)

###########GEO age

ttest1_str1 <- ttest1[,c(1,2,575)]
ttest2_str2 <- ttest2[,c(1,2,575)]
allltest_str <- rbind(ttest1_str1,ttest2_str2)

clinical30219 <- read.csv("clinical30219.csv")
clinical30219 <- clinical30219[,c("Geo_accession","Age.at.surgery")]
colnames(clinical30219) <- c("Geo_accession","Age")
clinical31210 <- read.csv("clinical31210.csv")
clinical31210 <- clinical31210[,c("Geo_accession","Age..years.")]
colnames(clinical31210) <- c("Geo_accession","Age")
clinical37735 <- read.csv("clinical37745.csv")
clinical37735 <- clinical37735[,c("Geo_accession","Age")]
clinical50081 <- read.csv("clinical50081.csv")
clinical50081 <- clinical50081[,c("Geo_accession","Age")]


clinicalallage <- rbind(clinical30219,clinical31210,clinical37735,clinical50081)
clinicalallage$age <- ifelse(clinicalallage$Age>65,'yes','no')
allltest_str$Geo_accession <- row.names(allltest_str)
geoallage<- merge(allltest_str,clinicalallage,by='Geo_accession')

MALE_SUR <- filter(geoallage,geoallage$age=='yes')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=98)","Low risk (n=111)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'GEO >65 patients (n=209)')
ggsave('GEO>65.pdf',width = 7,height = 7)

############STAGEI,II
clinical30219 <- read.csv("clinical30219.csv")
clinical30219 <- clinical30219[,c("Geo_accession","Stage")]
clinical31210 <- read.csv("clinical31210.csv")
clinical31210 <- clinical31210[,c("Geo_accession","Stage")]
clinical37735 <- read.csv("clinical37745.csv")
clinical37735 <- clinical37735[,c("Geo_accession","Stage")]
clinical50081 <- read.csv("clinical50081.csv")
clinical50081 <- clinical50081[,c("Geo_accession","Stage")]
clinicalallstage <- rbind(clinical30219,clinical31210,clinical37735,clinical50081)



allltest_str$Geo_accession <- row.names(allltest_str)
geoallstage<- merge(allltest_str,clinicalallstage ,by='Geo_accession')

MALE_SUR <- filter(geoallstage,geoallstage$Stage=='I')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=131)","Low risk (n=278)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'GEO stage I patients (n=409)')
ggsave('geoI.pdf',width = 7,height = 7)


#####################DCA
source("stdca.R") 
tcgaclinical <- read.csv("TCGAclinical.csv",check.names = F)
DCAtcga_riskscore <- ttrain[,c(1,2,575)]
DCAtcga_riskscore$X_PATIENT <- row.names(DCAtcga_riskscore)
tcgaclinical_DCA <- tcgaclinical[,c(1,11,12,15)]
geoallageDCA <- merge(DCAtcga_riskscore,tcgaclinical_DCA,by="X_PATIENT")
#DCA_tcgaall <- na.omit(DCA_tcgaall)
SrvDCA = Surv(geoallageDCA$times, geoallageDCA$status)
colnames(geoallageDCA)[5] <- 'age'
colnames(geoallageDCA)[6] <- 'gender'
colnames(geoallageDCA)[7] <- 'stage'
write.csv(geoallageDCA,'TCGArisknomo.csv')
geoallageDCA$age <- ifelse(geoallageDCA$age>65,'Yes','No')
geoallageDCA$age <- as.factor(geoallageDCA$age)
geoallageDCA$risk <- as.factor(geoallageDCA$risk)

coxmod1 = coxph(SrvDCA ~ age, data=geoallageDCA)
summary(coxmod1)
geoallageDCA$model1 = c(1 - (summary(survfit(coxmod1,newdata=geoallageDCA), times=1095)$surv))

coxmod2 = coxph(SrvDCA ~  gender, data=geoallageDCA)
summary(coxmod2)
geoallageDCA$model2 = c(1 - (summary(survfit(coxmod2,newdata=geoallageDCA), times=1095)$surv))

coxmod3 = coxph(SrvDCA ~ stage, data=geoallageDCA)
summary(coxmod3)
geoallageDCA$model3 = c(1 - (summary(survfit(coxmod3,newdata=geoallageDCA), times=1095)$surv))

coxmod4 = coxph(SrvDCA ~ stage+age+gender, data=geoallageDCA)
geoallageDCA$model4 = c(1 - (summary(survfit(coxmod4,newdata=geoallageDCA), times=1095)$surv))

coxmod5 = coxph(SrvDCA ~ risk, data=geoallageDCA)
geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
summary(coxmod5)
summary(coxmod5)
geoallageDCA$model5 = c(1 - (summary(survfit(coxmod5,newdata=geoallageDCA), times=1095)$surv))


##########
geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
geoallageDCA$age <- relevel(geoallageDCA$age,ref = 'No')
geoallageDCA$gender <- relevel(geoallageDCA$gender,ref = 'FEMALE')

coxmod6 = coxph(SrvDCA ~ risk+age+gender+stage, data=geoallageDCA)

summary(coxmod6)
geoallageDCA$model6 = c(1 - (summary(survfit(coxmod6,newdata=geoallageDCA), times=1095)$surv))


colnames(geoallageDCA)[8:13] <- c('Net Benefit: Age','Net Benefit: Gender','Net Benefit: Stage',
                                 'Net Benefit: Age+Gender+Stage','Net Benefit: Prediction model',
                                 'Net Benefit: Prediction model+Age+Gender+Stage')

pdf("net_benefit3-year_riskTCGA.pdf",width = 8,height = 8)
stdca(geoallageDCA, outcome="status", ttoutcome="times", timepoint=1095,
      predictors=c('Net Benefit: Age','Net Benefit: Gender','Net Benefit: Stage','Net Benefit: Age+Gender+Stage',
                   'Net Benefit: Prediction model',
                   'Net Benefit: Prediction model+Age+Gender+Stage'), smooth=TRUE,xstop = 0.5)
dev.off()


################HR TCGA
geoallageDCA$Age <- ifelse(geoallageDCA$age>65,'Yes','No')
geoallageDCA$Age <- as.factor(geoallageDCA$Age)
geoallageDCA$risk <- as.factor(geoallageDCA$risk)

geoallageDCA$Age <- relevel(geoallageDCA$Age,ref = 'No')
coxmod1 = coxph(SrvDCA ~ Age , data=geoallageDCA)
summary(coxmod1)
aaaa <- summary(coxmod1)
aaaa$coefficients[,5]
aaaa$conf.int[,c(1,3,4)]

geoallageDCA$gender <- relevel(geoallageDCA$gender,ref = 'FEMALE')
coxmod2 = coxph(SrvDCA ~  gender, data=geoallageDCA)
summary(coxmod2)

coxmod3 = coxph(SrvDCA ~ stage, data=geoallageDCA)
summary(coxmod3)

geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
coxmod5 = coxph(SrvDCA ~ risk, data=geoallageDCA)
summary(coxmod5)

########################TCGAnomogram
library(survival)
library(rms)
??datadist
geoallageDCA <- read.csv("georisknomo.csv",row.names = 1)
geoallageDCA <- geoallageDCA[,-1]
geoallageDCA$age <- ifelse(geoallageDCA$age>65,'Yes','No')
geoallageDCA$Stage1 <- str_replace_all(geoallageDCA$Stage,'II','I')
geoallageDCA$Stage2 <- str_replace_all(geoallageDCA$Stage1,'III','IV')
geoallageDCA <- geoallageDCA[,c(1,2,3,4,5,8)]


dd<-datadist(geoallageDCA)
dd
options(datadist="dd")
options(na.action="na.delete")
summary(geoallageDCA$times)

geoallageDCA <- merge(DCAtcga_riskscore,tcgaclinical_DCA,by="X_PATIENT")
geoallageDCA <- na.omit(geoallageDCA)
geoallageDCA$risk <- as.factor(geoallageDCA$risk)
geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
colnames(geoallageDCA)[4] <- 'Risk'
colnames(geoallageDCA)[5] <- 'Age'
colnames(geoallageDCA)[6] <- 'Gender'
colnames(geoallageDCA)[7] <- 'Stage'

############

coxpbc<-cph(formula = Surv(times,status) ~Risk+Age + Stage +Gender ,data=geoallageDCA,x=T,y=T,surv = T,na.action = na.delete)
print(coxpbc)
surv<-Survival(coxpbc) 

surv1<-function(x) surv(365,x)
surv3<-function(x) surv(1095,x)
surv5<-function(x) surv(1825,x)

x <- nomogram(coxpbc,fun = list(surv1,surv3,surv5),
            funlabel = c('1-Year OS Rate','3-Year OS Rate','5-Year OS Rate'),
            maxscale = 100)
pdf("nomogram_classical.pdf",width = 14,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.30,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.2, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = 'white')
dev.off()

f5<-cph(formula = Surv(times,status) ~Risk+Age + Stage +Gender ,data=geoallageDCA,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 
365*5
#参数m=50表示每组50个样本进行重复计算
cal5<-calibrate(f5, cmethod="KM", method="boot",u=1825,m=300,B=2000) 

mypal
pdf("calibration_5yall.pdf",width = 5,height = 5)
plot(cal5,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("red"),#error bar的颜色
     xlim = c(0.32,0.78),ylim= c(0.28,0.88),
     xlab = "Nomogram-prediced 5-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 1, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("gold")#对角线的颜色
) 
dev.off()







######################GEO
clinical30219 <- read.csv("clinical30219.csv")
clinical30219 <- clinical30219[,c(1,8,11,18)]
colnames(clinical30219) <- c('Geo_accession','age','gender','stage')

clinical31210 <- read.csv("clinical31210.csv")
clinical31210 <- clinical31210[,c(1,9,13,17)]
colnames(clinical31210) <- c('Geo_accession','age','gender','stage')


clinical37735 <- read.csv("clinical37745.csv")
clinical37735 <- clinical37735[,c(1,7,11,14)]
colnames(clinical37735) <- c('Geo_accession','age','gender','stage')
clinical50081 <- read.csv("clinical50081.csv")
clinical50081 <- clinical50081[,c(1,8,14,20)]
colnames(clinical50081 ) <- c('Geo_accession','age','gender','stage')

clinicalallstage <- rbind(clinical30219,clinical31210,clinical37735,clinical50081)
allltest_str
geoallageDCA <- merge(allltest_str,clinicalallstage,by='Geo_accession')
write.csv(geoallageDCA,'georisknomo.csv')

SrvDCA_geo= Surv(geoallageDCA$times, geoallageDCA$status)

coxmod1 = coxph(SrvDCA_geo ~ age , data=geoallageDCA)
summary(coxmod1)
geoallageDCA$model1 = c(1 - (summary(survfit(coxmod1,newdata=geoallageDCA), times=1025)$surv))


coxmod2 = coxph(SrvDCA_geo ~  gender, data=geoallageDCA)
geoallageDCA$model2 = c(1 - (summary(survfit(coxmod2,newdata=geoallageDCA), times=1025)$surv))

coxmod3 = coxph(SrvDCA_geo ~ stage, data=geoallageDCA)
geoallageDCA$model3 = c(1 - (summary(survfit(coxmod3,newdata=geoallageDCA), times=1025)$surv))

coxmod4 = coxph(SrvDCA_geo ~ stage+age+gender, data=geoallageDCA)
geoallageDCA$model4 = c(1 - (summary(survfit(coxmod4,newdata=geoallageDCA), times=1025)$surv))

geoallageDCA$risk <- as.factor(geoallageDCA$risk)
geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')

coxmod5 = coxph(SrvDCA_geo ~ risk, data=geoallageDCA)

summary(coxmod5)
geoallageDCA$model5 = c(1 - (summary(survfit(coxmod5,newdata=geoallageDCA), times=1025)$surv))

#########
geoallageDCA$risk <- as.factor(geoallageDCA$risk)
geoallageDCA$age <- as.factor(ifelse(geoallageDCA$age>65,'Yes','No'))

geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
geoallageDCA$status
  
coxmod6 = coxph(SrvDCA_geo ~ risk+age+gender+stage, data=geoallageDCA)
summary(coxmod6)



geoallageDCA$model6 = c(1 - (summary(survfit(coxmod6,newdata=geoallageDCA), times=1025)$surv))


colnames(geoallageDCA)[8:13] <- c('Net Benefit: Age','Net Benefit: Gender','Net Benefit: Stage',
                                  'Net Benefit: Age+Gender+Stage','Net Benefit: Prediction model',
                                  'Net Benefit: Prediction model+Age+Gender+Stage')

pdf("net_benefit3-year_riskGEO.pdf",width = 8,height = 8)
stdca(geoallageDCA, outcome="status", ttoutcome="times", timepoint=1025,
      predictors=c('Net Benefit: Age','Net Benefit: Gender','Net Benefit: Stage','Net Benefit: Age+Gender+Stage',
                   'Net Benefit: Prediction model+Age+Gender+Stage'), smooth=TRUE,xstop = 0.5)
dev.off()

####################GEO HR



geoallageDCA$Age <- ifelse(geoallageDCA$age>65,'Yes','No')
geoallageDCA$Age <- as.factor(geoallageDCA$Age)
geoallageDCA$risk <- as.factor(geoallageDCA$risk)

geoallageDCA$Age <- relevel(geoallageDCA$Age,ref = 'No')
coxmod1 = coxph(SrvDCA_geo ~ Age , data=geoallageDCA)
summary(coxmod1)


geoallageDCA$gender <- relevel(geoallageDCA$gender,ref = 'Female')
coxmod2 = coxph(SrvDCA_geo ~  gender, data=geoallageDCA)
summary(coxmod2)

coxmod3 = coxph(SrvDCA_geo ~ stage, data=geoallageDCA)
summary(coxmod3)

geoallageDCA$risk <- relevel(geoallageDCA$risk,ref = 'low')
coxmod5 = coxph(SrvDCA_geo ~ risk, data=geoallageDCA)
summary(coxmod5)

################################nomoall

library(glmnet)
x1 <- data.matrix(TCGALUAD_LASSO[,colnames(TCGALUAD_LASSO) %in% best.vars])
y1 <- data.matrix(Surv(TCGALUAD_LASSO[,]$times,TCGALUAD_LASSO[,]$status))
cv.fitA11 <- cv.glmnet(x1,y1,type.measure = "deviance",family='cox')
pdf('DVlasso.pdf',height = 7,width = 9)
plot(cv.fitA11)
dev.off()
cv.fitA11$

###########
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","red"),2)

#设置x轴最大值
xmax <- 1.75

plotCoef_plus <- function (beta, norm, lambda, df, dev, label = FALSE, legend = FALSE, xvar = c("norm", 
                                                                                                "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) 
{
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  
  if (legend){
    #在右侧留出画图例的地方
    par(xpd = T, mar = par()$mar + c(0,0,0,6))
  }
  
  #修改bty，换个更好看的边框，还可以改成，o / n / 7 / l / c / u / ]
  if (is.null(type)) 
    matplot(index, t(beta), lty = 1, lwd = 2,
            xlab = xlab, ylab = ylab, 
            xlim = c(0, xmax), #设置x轴最大值
            col = mycol,#线的颜色
            type = "l", cex.lab=1.2, cex.axis=1,
             ...)#不画右边框
  else matplot(index, t(beta), lty = 1, lwd = 2,
               xlab = xlab, ylab = ylab, 
               xlim = c(0, xmax), 
               col = mycol,
               type = "l", cex.lab=1.2, cex.axis=1,
               bty="n", ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    
    #原函数打印序号，修改为打印基因名
    text(xpos, ypos, paste(colnames(x1)[which]),
         cex = 0.5, #基因名字体大小
         #基因名的颜色跟线一样
         col = mycol,
         #如果你不想要彩色的字，就用下面这行
         #col = "black",
         pos = pos)
  }
  if (legend) {
    #画图例
    legend("topright",
           inset=c(-0.12,0),#图例画到图外面
           legend = rownames(myexpr), #图例文字
           col = mycol, #图例线的颜色，与文字对应
           lwd = 3, #图例中线的粗细
           cex = 1, #图例字体大小
           bty = "n") #不显示图例边框
  }
  par(xpd=FALSE)
}

plot.glmnet_plus <- function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, legend = FALSE,
                              ...) 
{
  xvar = match.arg(xvar)
  plotCoef_plus(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
                label = label, legend = legend, xvar = xvar, ...)
}
pdf("lasso.pdf",width = 9,height = 7)
plot.glmnet_plus(cv.fitA11$glmnet.fit, label = TRUE, #打印基因名
                 legend = FALSE) #不显示图例
dev.off()
best.coef
write.csv(best.vars,'bestbar.csv')
best.vars


################redo nomogram
library(rms)
tcga <- read.csv("tcganomo.csv",row.names = 1)
geo <- read.csv("geonomo.csv",row.names = 1)

dd<-datadist(tcga)
dd
options(datadist="dd")
coxpbc<-cph(formula = Surv(times,status) ~Risk+Age + Stage +Gender ,data=tcga,x=T,y=T,surv = T,na.action = na.delete)
print(coxpbc)
surv<-Survival(coxpbc) 
surv1<-function(x) surv(365,x)
surv3<-function(x) surv(1095,x)
surv5<-function(x) surv(1825,x)

x <- nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp = FALSE,
              funlabel = c('1-Year OS Rate','3-Year OS Rate','5-Year OS Rate'),
              maxscale = 100)
pdf('nomotcga.pdf',width = 13,height = 8)
plot(x)
dev.off()
##########################

allpatient <- rbind(tcga,geo)
allpatient$pre <- predict(coxpbc,newdata = allpatient)

##################

allcox<-cph(formula = Surv(times,status) ~pre ,data=allpatient,x=T,y=T,surv = T,na.action=na.delete,
            time.inc = 365)
allcox

Cindex <- rcorrcens(Surv(allpatient$times,allpatient$status)~predict(allcox))
Cindex

all<-calibrate(allcox, cmethod="KM", method="boot",u=365,m=300,B=20000)

v365 <- validate(allcox, method="boot",u=365,m=300,B=20000,dxy = TRUE)

Dxy365 = v365[rownames(v365)=='Dxy', colnames(v365)=='index.corrected']

bias_corrected_c_index365  <- abs(Dxy365)/2+0.5  # 计算校正c-index

bias_corrected_c_index365

pdf('all1year.pdf',height = 5,width = 5)
plot(all,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("red"),#error bar的颜色
     xlab = "Nomogram-prediced 1-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7) #字的大小
lines(all[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 1, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("gold")#对角线的颜色
) 
dev.off()

#################1095
allcox<-cph(formula = Surv(times,status) ~pre ,data=allpatient,x=T,y=T,surv = T,na.action=na.delete,
            time.inc = 1095)

all<-calibrate(allcox, cmethod="KM", method="boot",u=1095,m=300,B=20000)
v365 <- validate(allcox, method="boot",u=1095,m=300,B=20000,dxy = TRUE)
v365
Dxy365 = v365[rownames(v365)=='Dxy', colnames(v365)=='index.corrected']
bias_corrected_c_index1095  <- abs(Dxy365)/2+0.5  # 计算校正c-index

pdf('all5year.pdf',height = 5,width = 5)
plot(all,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("red"),#error bar的颜色
     xlab = "Nomogram-prediced 5-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7) #字的大小
lines(all[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 1, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("gold")#对角线的颜色
) 
dev.off()


##############1825
allcox<-cph(formula = Surv(times,status) ~pre ,data=allpatient,x=T,y=T,surv = T,na.action=na.delete,
            time.inc = 1825)
allcox

Cindex <- rcorrcens(Surv(allpatient$times,allpatient$status)~predict(allcox))
Cindex
all<-calibrate(allcox, cmethod="KM", method="boot",u=1825,m=300,B=20000)
v365 <- validate(allcox, method="boot",u=1825,m=300,B=20000,dxy = TRUE)
v365
Dxy365 = v365[rownames(v365)=='Dxy', colnames(v365)=='index.corrected']

bias_corrected_c_index1825 <- abs(Dxy365)/2+0.5  # 计算校正c-index

plot(all,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("red"),#error bar的颜色
     xlab = "Nomogram-prediced 5-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 1, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("gold")#对角线的颜色
) 




############tcga smoking



tcgaclinical <- read.csv("TCGAclinical.csv",check.names = F)

strtcga_stage <- ttrain[,c(1,2,575)]

tcgaclinical_smoking <- tcgaclinical[,c(1,24)]

strtcga_stage$X_PATIENT <- row.names(strtcga_stage)

str_STAGE_smoking <- merge(strtcga_stage,tcgaclinical_smoking,by="X_PATIENT")
str_STAGE_smoking$smoke <- ifelse(str_STAGE_smoking$tobacco_smoking_history==1,'Non-smoker','Ever-smoker')

smoking_SUR <- filter(str_STAGE_smoking,str_STAGE_smoking$smoke=='Non-smoker')

sum(smoking_SUR$risk=='high')
sum(smoking_SUR$risk=='low')

Sur_stage1 <- Surv(smoking_SUR$times/365,smoking_SUR$status)
sfit_stage1 <- survfit(Sur_stage1 ~ risk,data=smoking_SUR )
group <- factor(smoking_SUR $risk, levels = c("low", "high"))
data.survdiff <- survdiff(Sur_stage1 ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
ggsurvplot(sfit_stage1,
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=10)","Low risk (n=59)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD non-somker (n=69)')
ggsave('non-smokertcga.pdf',width = 7,height = 7)


#######GEO_smoking

clinical31210_somking <- clinical31210[,c("Geo_accession",'Smoking.status')]

clinical50081_smoking <- clinical50081[,c("Geo_accession","Smoking")]
colnames(clinical31210_somking)[2] <- 'Smoking'

clinical50081_smoking <- filter(clinical50081_smoking,clinical50081_smoking$Smoking!='Unable to determine')
clinical50081_smoking$Smoking <- ifelse(clinical50081_smoking$Smoking=='Never','Never-smoker','Ever-smoker')
clinicalallsmoking <- rbind(clinical31210_somking,clinical50081_smoking)


geoallstage
allltest_str$Geo_accession <- row.names(allltest_str)
geoallstage<- merge(allltest_str,clinicalallsmoking ,by='Geo_accession')

MALE_SUR <- filter(geoallstage,geoallstage$Smoking=='Ever-smoker')
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
           conf.int=F, #置信区间
           #fun="pct",
           legend = c(0.8,0.9),
           palette = c("#EFC000FF", "#0073C2FF"),
           xlab=('Years'),legend.title='',main=c('a'),
           pval.method = F,
           risk.table =F,
           ncensor.plot = F,
           legend.labs=c("High risk (n=84)","Low risk (n=119)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'GEO ever-smoker (n=203)')
ggsave('geom-eversmoker.pdf',width = 7,height = 7)


SrvDCA = Surv(str_STAGE_smoking$times,str_STAGE_smoking$status)
str_STAGE_smoking$smoke <- as.factor(str_STAGE_smoking$smoke)
str_STAGE_smoking$smoke <- relevel(str_STAGE_smoking$smoke,ref = 'Ever-smoker')

coxmod5 = coxph(SrvDCA ~smoke, data=str_STAGE_smoking)
summary(coxmod5)

geoallageDCA <- merge(geoallageDCA ,geoallstage[,c(1,5)],by='Geo_accession')

surr <- Surv(geoallageDCA$times.x,geoallageDCA$status.x)
geoallageDCA <- geoallageDCA[,c(1,2,3,4,5,6,7,11)]
coxxx <- coxph(surr~risk.x+age+gender+stage+Smoking.x,data = geoallageDCA)
summary(coxxx)
colnames(geoallageDCA )

