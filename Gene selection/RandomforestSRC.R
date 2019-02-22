library("randomForestSRC")
library("ggRandomForests")
library("survminer")
library("tidyverse")
library("plyr")
library("survival")
library("ggsci")
TCGALUAD<-read.csv('TCGALUAD.csv',row.name=1)
TCGALUAD <- TCGALUAD[TCGALUAD$times>=30,]     #follow-up days > 30
############################################
set.seed(2019)
TCGALUAD[1:4,1:4]
MultiNames<-colnames(TCGALUAD)[3:14479]
stp <- as.formula(paste0("Surv(times, status)~",paste(MultiNames,collapse = '+')))
stp
mtry<-sqrt(length(MultiNames)) 
nodesize<-3 # 3 for survival analysis
############################################
res.rsf <- rfsrc(stp,
                 data = TCGALUAD,
                 ntree = 1000,        
                 nodesize = nodesize,
                 splitrule = 'logrankscore',   #using logrank-score splitting method
                 forest = T,   
                 tree.err = T, 
                 importance = T,
                 mtry=mtry,
                 proximity=T, 
                 statistics = F)
############################################Variable selection 
vars <- var.select(object = res.rsf,
           cause =1,
           method = "vh.vimp", ########variable hunting using vh.vimp method
           conservative = c("medium"), # medium conservation 
           ntree = 1000,
           mvars = 3000, #randomly selecting 3000 genes for variable hunting
           mtry =  NULL,
           nodesize = nodesize, splitrule = "logrankscore", nsplit = 10, xvar.wt = NULL,
           refit = T, fast = T,
           na.action = c("na.impute"), 
           always.use = NULL, nrep = 100, K = 10, nstep = 200,  #repetiting 100 times
           prefit =  list(action = T, ntree = 1000,
                          mtry = mtry, nodesize = nodesize, nsplit = 1),
           verbose = TRUE)
vars$topvars.  #check top variable
