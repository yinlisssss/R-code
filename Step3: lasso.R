###write a function for splitting data into training, internal testing, and then validate in external cohorts.
cirlasso<-function(train=train,test1=test1,test2=test2,varlength= 30 ,nloop=100, modelcut=0.7,R2=0.7) {
  
  Train<-as.data.frame(train)
  Test1<-as.data.frame(test1)
  Test2<-as.data.frame(test2)
  commongene <- Reduce(intersect,list(colnames(Train)[3:ncol(Train)],colnames(Test1)[3:ncol(Test1)],colnames(Test2)[3:ncol(Test2)]))
  Train<-Train[,c("times","status", commongene)]
  Test1<-Test1[,c("times","status", commongene)]
  Test2<-Test1[,c("times","status", commongene)]
######
  FUN<-function(d,t1,t2,v,r){
    n<-nrow(d)
    res.tArc   <-c()
    res.tBrc   <-c()
    res.test1rc <-c()
    res.test2rc <-c()
    res.trainrc<-c()
    genelist   <-vector(mode="character",length=v)
    tArs       <-c()
    tBrs       <-c()
    test1rs     <-c()
    test2rs     <-c()
    trainrs    <-c()
    coef       <-vector(length=v)
    perm <- sample(1:n)
    tA <- perm[1:ceiling(r*n )]
    tB <- perm[(ceiling(r*n ) + 1):n]
    xA<-data.matrix(d[tA,3:length(colnames(d))])
    yA<-data.matrix(Surv(d[tA,]$times,d[tA,]$status))
    xB<-data.matrix(d[tB,3:length(colnames(d))])
    yB<-data.matrix(Surv(d[tB,]$times,d[tB,]$status))
    cv.fitA<-cv.glmnet(xA,yA,type.measure = "deviance",family = "cox")
    coeA<-coef(cv.fitA$glmnet.fit,s=cv.fitA$lambda.min)
    tA.active.coef<-coeA[which(coeA[,1]!=0)]
    tA.name<-row.names(coeA)[which(coeA[,1]!=0)]
    if(length(tA.name)>=3 & length(tA.name)<=v){
      tA.riskscore <- predict(cv.fitA$glmnet.fit, xA, s=cv.fitA$lambda.min, type="link")
      tArc<-rcorr.cens(-tA.riskscore, 
                       Surv(d[tA,]$times,d[tA,]$status))["C Index"] 
      tB.riskscore <- predict(cv.fitA$glmnet.fit, xB, s=cv.fitA$lambda.min, type="link")
      tBrc<-rcorr.cens(-tB.riskscore, 
                       Surv(d[tB,]$times,d[tB,]$status))["C Index"] 
      
      test1.riskscore <- predict(cv.fitA$glmnet.fit, as.matrix(t1[,3:ncol(t1)]), s=cv.fitA$lambda.min, type="link")
      test1.rc<-rcorr.cens(-test1.riskscore, 
                           Surv(t1$times,t1$status))["C Index"]
      
      test2.riskscore <- predict(cv.fitA$glmnet.fit, as.matrix(t2[,3:ncol(t2)]), s=cv.fitA$lambda.min, type="link")
      test2.rc<-rcorr.cens(-test2.riskscore, 
                           Surv(t2$times,t2$status))["C Index"]
      
      
      train.riskscore <- predict(cv.fitA$glmnet.fit, as.matrix(d[,3:ncol(d)]), s=cv.fitA$lambda.min, type="link")
      train.rc<-rcorr.cens(-train.riskscore, 
                           Surv(d$times,d$status))["C Index"]
      res.tArc<-rbind(res.tArc, tArc) #internal train c-index
      res.tBrc<-rbind(res.tBrc, tBrc) #internal test c-index
      res.test1rc<-rbind(res.test1rc, test1.rc) #external test c-index
      res.test2rc<-rbind(res.test2rc, test2.rc)
      res.trainrc<-rbind(res.trainrc,train.rc) #whole datasets c-index
      genelist<-rbind(genelist,tA.name) #gene
      tArs<-cbind(tArs, tA.riskscore)   #risk score for  train
      tBrs<-cbind(tBrs, tB.riskscore)   #risk score for internal test
      test1rs<-cbind(test1rs,test1.riskscore) #risk score for external test
      test2rs<-cbind(test2rs,test2.riskscore) #risk score for external validation 
      
      coef<-cbind(coef,tA.active.coef) 
      
      trainrs<-cbind(trainrs,train.riskscore) #all datasets risk score
    }
    cirlasso_results<-list(res.tArc=res.tArc,
                           res.tBrc=res.tBrc,
                           res.test1rc=res.test1rc,
                           res.test2rc=res.test2rc,
                           res.trainrc=res.trainrc,
                           genelist=genelist,
                           tArs=tArs,
                           tBrs=tBrs,
                           test1rs=test1rs,
                           test2rs=test2rs,
                           coef=coef,
                           trainrs=trainrs)
  }
  #parallel
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  clusterExport(cl,c("Surv","rcorr.cens","cv.glmnet","coef","predict"))
  mm <- foreach(i=1:nloop, .combine=rbind,.multicombine=F,.verbose=T )%dopar% {
    FUN(d=Train,t1=Test1,t2 = Test2,v=varlength,r=R2)
  }
  stopCluster(cl)
  #extract results
  options(warn =1)
  if(nloop==1){
    if(is.null(mm)){
      stop(paste("The nloops with ", nloop, " are too small!", sep=""))
    }else{
      res.tArc<-unlist(mm[1])
      res.tBrc<-unlist(mm[2])
      res.test1rc<-unlist(mm[3])
      res.test2rc<-unlist(mm[4])
      res.trainrc<-unlist(mm[5])
      genelist<-mm[[6]]
      genelist<-genelist[substr(rownames(genelist),1,7)=="tA.name",]
      genelist<-t(as.matrix(genelist))
      tArs<-mm[[7]]
      tBrs <-mm[[8]]
      test1rs<-mm[[9]]
      test2rs<-mm[[10]]
      coef <-mm[[11]]
      coef<-coef[,!colnames(coef)=="coef"]
      coef<-as.matrix(coef)
      trainrs <-mm[[12]]
    }
  }else{
    if(is.null(mm)){
      stop(paste("The nloops with ", nloop, " are too small!", sep=""))
    }else{
      res.tArc<-do.call("rbind", mm[,1])
      res.tBrc<-do.call("rbind", mm[,2])
      res.test1rc<-do.call("rbind", mm[,3])
      res.test2rc<-do.call("rbind", mm[,4])
      res.trainrc<-do.call("rbind", mm[,5])
      genelist<-do.call("rbind", mm[,6])
      genelist<-genelist[substr(rownames(genelist),1,7)=="tA.name",]
      if(is.null(dim(genelist))){
        genelist<-t(genelist) 
      }else{
        genelist<-genelist
      }
      tArs<-do.call("cbind", mm[,7])
      tBrs<-do.call("cbind", mm[,8])
      test1rs<-do.call("cbind", mm[,9])
      test2rs<-do.call("cbind", mm[,10])
      coef<-do.call("cbind", mm[,11])
      coef<-as.matrix(coef[,colnames(coef)=="tA.active.coef"])
      trainrs<- do.call("cbind", mm[,12])
    }
  }
  

  rcindex<-data.frame(res.tArc,res.tBrc,res.test1rc,res.test2rc,res.trainrc,row.names = 1:length(res.tArc))
  colnames(rcindex)<-c("res.tArc.cindex","res.tBrc.cindex","res.test1rc.cindex","res.test2rc.cindex","res.trainrc.cindex")
  rcindex$index<-rcindex$res.test1rc.cindex^2+rcindex$res.test2rc.cindex^2+rcindex$res.tBrc.cindex^2 #变量筛选指数index的定义
  colnames(tArs)<- 1:ncol(tArs)
  colnames(tBrs)<- 1:ncol(tBrs)
  colnames(test1rs)<- 1:ncol(test1rs)
  colnames(test2rs)<- 1:ncol(test2rs)
  colnames(trainrs)<- 1:ncol(trainrs)
  colnames(coef)<- 1:ncol(coef)
  rownames(genelist)<-c(1:nrow(genelist))
  

  selct.rcindex<-rcindex[rcindex$res.trainrc.cindex>modelcut & rcindex$res.test1rc.cindex>modelcut & rcindex$res.test2rc.cindex>modelcut,]#tA和tB，test里c-index均大于0.7
  options(warn =1)
  if (is.na(selct.rcindex[1,1])) { 
    warning(paste0("The modelcut with ", modelcut, " are too high!"," Intermediate results are saved in cirlasso_intermediate_results.rda "))
    save(rcindex,
         selct.rcindex,
         genelist,
         tArs,
         tBrs,
         test1rs,
         test2rs,
         trainrs,
         coef,
         file = "cirlasso_intermediate_results.rda") 
    
    cirlasso_intermediate_results<-list(rcindex=rcindex,
                                        genelist=genelist,
                                        tArs=tArs,
                                        tBrs=tBrs,
                                        test1rs=test1rs,
                                        test2rs=test2rs,
                                        trainrs=trainrs,
                                        coef=coef)
    
    
    
  }else{
    best.rcindex<-selct.rcindex[selct.rcindex$index==max(selct.rcindex$index),] 
    best.vars<-unique(genelist[rownames(best.rcindex),]) 
    best.tArs<-tArs[,rownames(best.rcindex)]
    names(best.tArs)<-rownames(tArs)
    best.tBrs<-tBrs[,rownames(best.rcindex)]
    names(best.tBrs)<-rownames(tBrs)
    best.test1rs<-test1rs[,rownames(best.rcindex)]
    names(best.test1rs)<-rownames(test1rs)
    best.test2rs<-test2rs[,rownames(best.rcindex)]
    names(best.test2rs)<-rownames(test2rs)
    best.coef<-unique(coef[,rownames(best.rcindex)])
    names(best.coef)<-c(best.vars)
    best.trrainrs<-trainrs[,rownames(best.rcindex)]
    names(best.trrainrs)<-rownames(trainrs)
    save(rcindex,
         genelist,
         tArs,
         tBrs,
         test1rs,
         test2rs,
         trainrs,
         best.rcindex,
         best.vars,
         best.tArs,
         best.tBrs,
         best.test1rs,
         best.test2rs,
         best.coef,
         best.trrainrs,
         coef,
         file = "cirlasso_results.rda") #保存结果
    cirlasso_results<-list(rcindex=rcindex,
                           genelist=genelist,
                           tArs=tArs,
                           tBrs=tBrs,
                           test1rs=test1rs,
                           test2rs=test2rs,
                           trainrs=trainrs,
                           best.rcindex=best.rcindex,
                           best.vars=best.vars,
                           best.tArs=best.tArs,
                           best.tBrs=best.tBrs,
                           best.test1rs=best.test1rs,
                           best.test2rs=best.test2rs,
                           best.coef=best.coef,
                           best.trrainrs=best.trrainrs,
                           coef=coef)
    return(cirlasso_results)
  }
  
}

###################varlength: variables should be less than 40,
###################nloop: number of loop, 100,000 times
###################modelcut: model auc > 0.68
###################R2=0.7: randomly 70% to training, 30% to internal testing sets
laresults<-cirlasso(train=TCGALUAD,
                    test1=external_test,
                    test2=external_validation,varlength=40,nloop=100000, modelcut=0.68, R2=0.7)





