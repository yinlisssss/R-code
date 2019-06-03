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
