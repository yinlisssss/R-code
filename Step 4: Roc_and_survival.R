
########Once Lasso 100,000 loops done, draw roc curve and survival analysis

ttrain<-TCGALUAD_LASSO 
ttrain$riskscore <- best.trrainrs
res.cutsur <- surv_cutpoint(ttrain, time = "times",  #use best separation (survminer package)
                         event = "status", 
                         variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]
risk<-as.vector(ifelse(best.trrainrs>=res.cut,"high","low"))
ttrain$risk<-risk
sum(ttrain$risk=="high")
sum(ttrain$risk=="low")
write.csv(ttrain[,c(574,575)],'tcgarisk.csv')

# TCGA set
Sur <- Surv(ttrain$times/365,ttrain$status)
sfit <- survfit(Sur ~ risk,data=ttrain)
group <- factor(ttrain$risk, levels = c("low", "high"))
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
           legend.labs=c("High risk (n=79)","Low risk (n=413)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                           paste("p = ",round(p.val,3), sep = "")),
             HR, CI, sep = "\n"))+labs(title = 'TCGA LUAD cohort (n=492)')
ggsave("tcga.pdf",width = 7,height = 7)


# ROC 

survivalROC_helper <- function(t,data) {
  survivalROC(Stime=data$times/365, status=data$status, marker = data$riskscore,
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  mutate(survivalROC = map(t, data=ttrain, survivalROC_helper),
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
  scale_color_jco(labels=c("AUC at 1-year: 0.753","AUC at 3-year: 0.726","AUC at 5-year: 0.656"))+
  guides(color=guide_legend(title=NULL))
ggsave('trainroc.pdf',height = 7,width = 7)


#c factor plot
fp<- best.trrainrs
names(fp)<-rownames(trainrs)
fp_dat<-data.frame(s=1:length(fp),v=as.numeric(sort(fp)))
fp_dat$Risk<-ifelse(fp_dat$v>=res.cut,"high","low")
sur_dat<-data.frame(s=1:length(fp),
                    t=ttrain[names(sort(fp)),'times']/365,
                    e=ttrain[names(sort(fp)),'status'])
sur_dat$Status<-as.factor(ifelse(sur_dat$e==0,'alive','death'))
exp_dat<-ttrain[names(sort(fp)),which(colnames(ttrain) %in% unique(as.matrix(best.vars)))]
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
  labs(y="Risk score",fill="Risk",title ='TCGA LUAD cohort (n=492)' )+
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
plot.h <- pheatmap(tmp,show_colnames = F,cluster_cols = F,
         color = colorRampPalette(c("#0073C2FF", "white", "#EFC000FF"))(100))

A <- plot_grid(plot.point, plot.sur,plot.h$gtable,
          labels = c("A", "",""),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,axis="t")
A
ggsave('trainrisk.pdf',height = 10,width = 8)




#8.2 external test
ttest1<-test1_lasso
ttest1$riskscore <- best.test1rs
res.cutsur <- surv_cutpoint(ttest1, time = "times", 
                            event = "status", 
                            variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]

ttest1.risk<-as.vector(ifelse(best.test1rs>=res.cut,"high","low"))
ttest1$risk<-ttest1.risk
sum(ttest1$risk=="high")
sum(ttest1$risk=="low")
write.csv(ttest1[,c(574,575)],'geo1risk.csv')
#a 
Sur <- Surv(ttest1$times/365,ttest1$status)
sfit <- survfit(Sur ~ ttest1.risk,data=ttest1)

group <- factor(ttest1$risk, levels = c("low", "high"))
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
           legend.labs=c("High risk (n=110)","Low risk (n=122)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'External testing set (n=232)')
ggsave('test1sur.pdf',height = 7,width = 7)




fp<- best.test1rs
names(fp)<-rownames(test1rs)
fp_dat<-data.frame(s=1:length(fp),v=as.numeric(sort(fp)))
fp_dat$Risk<-ifelse(fp_dat$v>=res.cut,"high","low")
sur_dat<-data.frame(s=1:length(fp),
                    t=ttest1[names(sort(fp)),'times']/365,
                    e=ttest1[names(sort(fp)),'status'])
sur_dat$Status<-as.factor(ifelse(sur_dat$e==0,'alive','death'))
exp_dat<-ttest1[names(sort(fp)),which(colnames(ttest1) %in% unique(as.matrix(best.vars)))]
plot.point<-ggplot(fp_dat,aes(x=s,y=v))+
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
  labs(y="Risk score",fill="Risk",title ='External testing set (n=232)' )+
  scale_color_manual(values = c("#EFC000FF", "#0073C2FF"),
                     labels=c("High","Low"))+guides(color=guide_legend(title=NULL))
plot.point

plot.sur<-ggplot(sur_dat,aes(x=s,y=t))+
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
plot.h<-pheatmap(tmp,show_colnames = F,cluster_cols = F,
                 color = colorRampPalette(c("#0073C2FF", "white", "#EFC000FF"))(100))
B <- plot_grid(plot.point, plot.sur,plot.h$gtable,
          labels = c("B", "",""),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,axis="b")
B

#external validation

ttest2<-test2_lasso
ttest2$riskscore <- best.test2rs
res.cutsur <- surv_cutpoint(ttest2, time = "times", 
                            event = "status", 
                            variables = 'riskscore')
res.cut <- res.cutsur$cutpoint[[1]]
ttest2.risk<-as.vector(ifelse(best.test2rs>=res.cut,"high","low"))
ttest2$risk<-ttest2.risk
sum(ttest2$risk=="high")
sum(ttest2$risk=="low")
write.csv(ttest2[,c(574,575)],'geo2risk.csv')

Sur <- Surv(ttest2$times/365,ttest2$status)
sfit <- survfit(Sur ~ ttest2.risk,data=ttest2)

group <- factor(ttest2$risk, levels = c("low", "high"))
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
           legend.labs=c("High risk (n=110)","Low risk (n=122)"),
           pval = paste(pval = ifelse(p.val < 0.0001, "p < 0.0001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))+labs(title = 'External testing set (n=232)')






fp<- best.test2rs
names(fp)<-rownames(test2rs)
fp_dat<-data.frame(s=1:length(fp),v=as.numeric(sort(fp)))
fp_dat$Risk<-ifelse(fp_dat$v>=res.cut,"high","low")
sur_dat<-data.frame(s=1:length(fp),
                    t=ttest2[names(sort(fp)),'times']/365,
                    e=ttest2[names(sort(fp)),'status'])
sur_dat$Status<-as.factor(ifelse(sur_dat$e==0,'alive','death'))
exp_dat<-ttest2[names(sort(fp)),which(colnames(ttest2) %in% unique(as.matrix(best.vars)))]
plot.point<-ggplot(fp_dat,aes(x=s,y=v))+
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
  labs(y="Risk score",fill="Risk",title ='External validation set (n=347)' )+
  scale_color_manual(values = c("#EFC000FF", "#0073C2FF"),
                     labels=c("High","Low"))+guides(color=guide_legend(title=NULL))
plot.point

plot.sur<-ggplot(sur_dat,aes(x=s,y=t))+
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
plot.h<-pheatmap(tmp,show_colnames = F,cluster_cols = F,
                 color = colorRampPalette(c("#0073C2FF", "white", "#EFC000FF"))(100))
C <- plot_grid(plot.point, plot.sur,plot.h$gtable,
               labels = c("C", "",""),
               label_x=0,
               label_y=1,
               align = 'v',ncol = 1,axis="b")
plot_grid(A,B,C,ncol = 3)

ggsave('riskall.pdf',height = 10,width = 20)








