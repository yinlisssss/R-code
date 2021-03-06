######Data acquisition and preprocessing
###TCGA
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
###clinical information: (https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018)
query.exp.down <- GDCquery(project = "TCGA-LUAD",
                           legacy = TRUE,
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq",
                           sample.type = c("Primary solid Tumor", "Solid Tissue Normal"),
                           experimental.strategy = "RNA-Seq")
GDCdownload(query.exp.down)
LUAD.exp <- GDCprepare(query = query.exp.down,
                       save = TRUE,
                       save.filename = "LUAD_RNAseq.rda")
dataPrep <- TCGAanalyze_Preprocessing(object = LUAD.exp, cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.25)


#Differential analysis, check all genes' log fold change
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
typesample = c("TP"))
# Diff.expr.analysis (DEA)
tcga <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 1 ,
                            logFC.cut = 0,
                            method = "glmLRT")
######GEO 
#GEO CEL files can be manually downloaded from https://www.ncbi.nlm.nih.gov/geo/ 
library(affy)
Data <- ReadAffy() ##read data in working directory
eset <- rma(Data)
write.exprs(eset, file="GES####.csv")#file name
###############
#Once all GEO data have been downloaded, merging all to remove batch effect.
library(sva)
library(limma)

##clinical information
library(GEOquery)
library(tidyverse)
library(rlang)
library(impute)
library(Biobase)
#for example use 31210, use GEOquery to download clinical information
gset <- getGEO("GSE31210", GSEMatrix =TRUE, getGPL = TRUE, AnnotGPL = TRUE,destdir = './')
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
exprdf<-data.frame(Biobase::exprs(gset))
dim(exprdf)
clinicaldata_31210 <- pData(gset)

#then you can extract the clinical information for subsequent analysis

#####read local files
GSE19188<- read.csv("GSE19188cel.csv",row.names = 1,header = T)
GSE30219<- read.csv("GSE30219cel.csv",row.names = 1,header = T)
GSE31210<- read.csv("GSE31210cel.csv",row.names = 1,header = T)
GSE37745<- read.csv("GSE37745cel.csv",row.names = 1,header = T)
GSE50081<- read.csv("GSE50081cel.csv",row.names = 1,header = T)
##
five_batch_cel <- Reduce(function(x,y) merge(x=x,y=y,by='gene'),list(GSE19188,GSE30219,GSE31210,GSE37745,GSE50081))
####sample characteristics
batchType_five=c(rep(1,110),rep(2,99),rep(3,246),rep(4,105),rep(5,127))
modType_five=c(rep("normal",65),rep("tumor",45),rep("normal",14),rep("tumor",85),
               rep("normal",20),rep("tumor",226),rep("normal",0),rep("tumor",105),
               rep("normal",0),rep("tumor",127))
mod_five = model.matrix(~as.factor(modType_five))
outTab_five=ComBat(five_batch_cel, batchType_five, mod_five, par.prior=TRUE)
finaldf_five <- as.data.frame(outTab_five)
#check distribution
boxplot(finaldf_five  ,col = "blue",xaxt = "n",outline = F,main="After batch normalization")
###normalization
rt1_five <- normalizeBetweenArrays(as.matrix(finaldf_five))
########
##data have been downloaded and pre-processed
