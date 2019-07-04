######use ICGC LUAD dataset to demonstrate that there is not overfitting.
####Download data from UCSC Xena
####https://xenabrowser.net/datapages/?cohort=ICGC%20(donor%20centric)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443


library(tidyverse)
library(data.table)
library(ggsci)
library(survminer)
library(survival)
library(survivalROC)

options(stringsAsFactors = FALSE)
###read local data
icgc_surival <- fread('donor.all_projects.survival.xena')
icgc_type <- fread('donor.all_projects.xena')

icgc_luad <- filter(icgc_type,icgc_type$project_code %like% 'LUAD')

icgc_luad <- icgc_luad[,c(1,2,3,4,9,10,11)]
icgc_surival_luad <- filter(icgc_surival,icgc_surival$icgc_donor_id %in% icgc_luad$sampleID)

icgc_surival_luad <- na.omit(icgc_surival_luad)
icgc_surival_luad <- filter(icgc_surival_luad,icgc_surival_luad$`_TIME_TO_EVENT`>=30)

##read lasso index (Supplementary table 5)
lasso_index <- fread('lasso.txt') 

lasso_index <- as.data.frame(lasso_index)
row.names(lasso_index) <- lasso_index[,1]
###read expression data
icgc_exp <- fread('exp_seq_donor_US_log2')
icgc_exp <- as.data.frame(icgc_exp)
row.names(icgc_exp) <- icgc_exp$donor
icgc_exp <- icgc_exp[,-1]
icgc_exp[1:4,1:4]

icgc_exp <- icgc_exp[lasso_index$Gene,]

aaa <- intersect(icgc_surival_luad$icgc_donor_id,colnames(icgc_exp))
icgc_exp_luad <- icgc_exp[,aaa]
icgc_exp_luad <- as.data.frame(t(icgc_exp_luad))
icgc_exp_luad <- as.matrix(icgc_exp_luad)
######extract LUAD clinical information
lasso_index <- lasso_index[,-1]

riskscore <- icgc_exp_luad %*% lasso_index 

riskscore <- as.data.frame(riskscore)
colnames(riskscore) <- 'riskscore'

riskscore$sampleid <- row.names(riskscore)

icgc_surival_luad$sampleid <- icgc_surival_luad$icgc_donor_id

icgc_risk_luad <- merge(icgc_surival_luad,riskscore,by='sampleid')


icgc_risk_luad <- icgc_risk_luad[,c(1,3,4,6)]

icgc_risk_luad$time <- icgc_risk_luad$`_TIME_TO_EVENT`/365
icgc_risk_luad$status <- ifelse(icgc_risk_luad$`_EVENT`=='alive',0,1)

icgc_risk_luad <- icgc_risk_luad[,c(1,5,6,4)]



######

###find the optimal cutpoint
res.cutsur <- surv_cutpoint(icgc_risk_luad , time = "time", 
                            event = "status", 
                            variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]


icgc_risk_luad$risk<-as.vector(ifelse(icgc_risk_luad$riskscore>=res.cut,"high","low"))

sum(icgc_risk_luad$risk=="high")
sum(icgc_risk_luad$risk=="low")

######surivavl analysis

Sur <- Surv(icgc_risk_luad$time,icgc_risk_luad$status)
sfit <- survfit(Sur ~ risk,data=icgc_risk_luad)
group <- factor(icgc_risk_luad$risk, levels = c("low", "high"))
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
           legend.labs=c("High risk (n=82)","Low risk (n=291)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'ICGC LUAD cohort (n=373)')
ggsave("icgc.pdf",width = 7,height = 7)


#####ROC curve

survivalROC_helper <- function(t,data) {
  survivalROC(Stime=data$time, status=data$status, marker = data$riskscore,
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  mutate(survivalROC = map(t, data=icgc_risk_luad, survivalROC_helper),
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
survivalROC_data1 <- survivalROC_data %>%
  mutate(auc =sprintf("%.3f",auc))%>%
  unite(month, t,auc,sep = " - year AUC: ")
AUC <-factor(survivalROC_data1$month)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2),element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(x = "1-Specificity",y="Sensitivity",title = 'Time-dependent ROC')+
  guides(color=guide_legend(title=NULL))
ggsave('icgcroc.pdf',height = 7,width = 7)



row.names( icgc_risk_luad) <- icgc_risk_luad$sampleid
fp<- icgc_risk_luad$riskscore
names(fp)<-rownames(icgc_risk_luad)
fp_dat<-data.frame(s=1:length(fp),v=as.numeric(sort(fp)))
fp_dat$Risk<-ifelse(fp_dat$v>=res.cut,"high","low")
sur_dat<-data.frame(s=1:length(fp),
                    t=icgc_risk_luad[names(sort(fp)),'time'],
                    e=icgc_risk_luad[names(sort(fp)),'status'])
sur_dat$Status<-as.factor(ifelse(sur_dat$e==0,'alive','death'))

exp_dat<-icgc_exp_luad[names(sort(fp)),]

plot.point <- ggplot(fp_dat,aes(x=s,y=v))+
  geom_point(aes(col=Risk),size=1)+
  geom_segment(aes(x = sum(fp_dat$Risk=="low"),
                   y = min(fp_dat$v),
                   xend = sum(fp_dat$Risk=="low"),
                   yend = res.cut),linetype="dashed")+
  geom_segment(aes(x=0,y=res.cut,
                   xend=sum(fp_dat$Risk=="low"),
                   yend=res.cut),linetype="dashed")+
  geom_text(aes(x=sum(fp_dat$Risk=="low")/1.5,
                y=res.cut+0.1,
                label=paste0("Cut-off value: ",round(res.cut,3))),size = 4,color="red")+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  labs(y="Risk score",fill="Risk",title ='ICGC LUAD cohort (n=373)' )+
  scale_color_manual(values = c("#EFC000FF", "#0073C2FF"),
                     labels=c("High","Low"))+guides(color=guide_legend(title=NULL))
plot.point

plot.sur <- ggplot(sur_dat,aes(x=s,y=t))+
  geom_point(aes(col=Status),size=1)+
  geom_vline(aes(xintercept=sum(fp_dat$Risk=="low")),linetype="dashed")+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_jco(labels=c("Alive", "Dead"))+
  labs(y="Survival time (Years)")+guides(color=guide_legend(title=NULL))
plot.sur 

tmp<-t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
plot.h <- pheatmap(tmp,show_colnames = F,cluster_cols = F,scale = 'row',
                   color = colorRampPalette(c("#0073C2FF", "white", "#EFC000FF"))(100))

A <- plot_grid(plot.point, plot.sur,plot.h$gtable,
               labels = c("C", "",""),
               label_x=0,
               label_y=1,
               align = 'v',ncol = 1,axis="t")
A

ggsave('icgcrisk.pdf',height = 10,width = 8)



#Stratification analysis
###########male

geneder_sample <- icgc_luad[,c(1,2,4)]
geneder_sample$sampleid <- geneder_sample$sampleID 
geneder_sample$age <- ifelse(geneder_sample$`_age_at_diagnosis`>65,'Yes','No')

male <- filter(geneder_sample,geneder_sample$`_gender`=='male')


malesur <- merge(icgc_risk_luad,male,by='sampleid')

sum(malesur $risk=="high")
sum(malesur $risk=="low")

Sur <- Surv(malesur$time,malesur$status)
sfit <- survfit(Sur ~ risk,data=malesur)
group <- factor(malesur$risk, levels = c("low", "high"))

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
           legend.labs=c("High risk (n=43)","Low risk (n=129)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'ICGC LUAD male cohort (n=172)')

ggsave("icgcmale.pdf",width = 7,height = 7)


###female
female <- filter(geneder_sample,geneder_sample$`_gender`=='female')

femalesur <- merge(icgc_risk_luad,female,by='sampleid')

sum(femalesur $risk=="high")
sum(femalesur $risk=="low")

Sur <- Surv(femalesur$time,femalesur$status)
sfit <- survfit(Sur ~ risk,data=femalesur)
group <- factor(femalesur$risk, levels = c("low", "high"))

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
           legend.labs=c("High risk (n=39)","Low risk (n=162)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'ICGC LUAD female cohort (n=201)')

ggsave("icgcfemale.pdf",width = 7,height = 7)

##########age65



age65 <- filter(geneder_sample,geneder_sample$age=='Yes')

age65 <- merge(icgc_risk_luad,age65,by='sampleid')

sum(age65 $risk=="high")
sum(age65 $risk=="low")

Sur <- Surv(age65$time,age65$status)
sfit <- survfit(Sur ~ risk,data=age65)
group <- factor(age65$risk, levels = c("low", "high"))

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
           legend.labs=c("High risk (n=47)","Low risk (n=161)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'ICGC LUAD >65 patients (n=208)')

ggsave("icgcage65.pdf",width = 7,height = 7)

####age<65
age65 <- filter(geneder_sample,geneder_sample$age=='No')

age65 <- merge(icgc_risk_luad,age65,by='sampleid')

sum(age65 $risk=="high")
sum(age65 $risk=="low")

Sur <- Surv(age65$time,age65$status)
sfit <- survfit(Sur ~ risk,data=age65)
group <- factor(age65$risk, levels = c("low", "high"))

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
           legend.labs=c("High risk (n=35)","Low risk (n=130)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'ICGC LUAD <=65 patients (n=165)')

ggsave("icgcageless65.pdf",width = 7,height = 7)








