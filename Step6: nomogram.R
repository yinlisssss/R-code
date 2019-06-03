tcga <- read.csv("tcganomo.csv",row.names = 1)
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

###########calibration
geo <- read.csv("geonomo.csv",row.names = 1)
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
     lwd = 2,
     lty = 1,
     errbar.col = c("red"),
     xlab = "Nomogram-prediced 1-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
lines(all[,c('mean.predicted',"KM")], 
      type = 'b',
      lwd = 2, 
      pch = 16, 
      col = c("red")) 
mtext("")
box(lwd = 1) 
abline(0,1,lty = 1,
       lwd = 2,
       col = c("gold")
) 
dev.off()

#################1095
allcox<-cph(formula = Surv(times,status) ~pre ,data=allpatient,x=T,y=T,surv = T,na.action=na.delete,
            time.inc = 1095)

all<-calibrate(allcox, cmethod="KM", method="boot",u=1095,m=300,B=20000)
v365 <- validate(allcox, method="boot",u=1095,m=300,B=20000,dxy = TRUE)
v365
Dxy365 = v365[rownames(v365)=='Dxy', colnames(v365)=='index.corrected']
bias_corrected_c_index1095  <- abs(Dxy365)/2+0.5  

pdf('all5year.pdf',height = 5,width = 5)
plot(all,
     lwd = 2,
     lty = 1,
     errbar.col = c("red"),
     xlab = "Nomogram-prediced 5-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
lines(all[,c('mean.predicted',"KM")], 
      type = 'b', 
      lwd = 2, 
      pch = 16,
      col = c("red")) 
mtext("")
box(lwd = 1) 
abline(0,1,lty = 1,
       lwd = 2,
       col = c("gold")
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
     lwd = 2,
     lty = 1,
     errbar.col = c("red"),
     xlab = "Nomogram-prediced 5-Year OS (%) In All Patients",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) 
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b', 
      lwd = 2,
      pch = 16,
      col = c("red")) 
mtext("")
box(lwd = 1) 
abline(0,1,lty = 1, 
       lwd = 2, 
       col = c("gold")
) 



