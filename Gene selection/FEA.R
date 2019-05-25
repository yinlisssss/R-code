#Functional enrichment analysis using moonlightR

devtools::install_github(repo = "ibsquare/MoonlightR")
library(MoonlightR)

##tcga can be obtained in step 1 data_dowload.R 

deg <- tcga[Genes_identified_by_three_methods,]

dataFEA <- FEA(DEGsmatrix = deg,
               BPname = names(DiseaseList))
dataFEA_filt <- dataFEA[dataFEA$FDR < 0.05,]
dataFEA_filt <- dataFEA_filt[abs(dataFEA_filt$Moonlight.Z.score) > 1,]


plotFEA1(dataFEA = dataFEA_filt,
        topBP = nrow(dataFEA_filt),
        height = 10,
        width = 15,
        additionalFilename = "TCGA_LUAD.pdf")


write.csv(dataFEA_filt,'FEA.csv')

