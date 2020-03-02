rm(list=ls())
setwd("/home/coen/faerslincs/")

library(caret)
library(glmnet)
library(doParallel)
library(dplyr)
library(reshape2)

source("functions.R")

# load genes original signature
load("results/ALLRES_2018-08-24__regorafenib-sunitinib_nfold_5_nbs_1000_ndesc_26_.Rdata")
signature<- vip1$gene


# FAERS new
# FAERS risk scores 
load("FAERS_external/riskscoresMetaJuly.Rdata") # generated in /data/lincs_predictCT/estimateCTrisk_finalJuly11.R
res<-resFAERS
res$drug<-as.character(res$drug)
names(res)[2]<-"zscore" # rename to zscore because some functions expect this naming
res<-select(res,drug,zscore,CIup,CIdown)
res$ae_descr="CT"
reflist2<-res
reflist_ext<-reflist2[reflist2$drug %in%  c("ibrutinib",  "lenvatinib",   "nintedanib"),]



# old risk scores
# FAERS risk scores 
load("FAERS_external/riskscoresMetaOrig.Rdata")
res<-resMeta
res$drug<-as.character(res$drug)
res<-select(res,drug,rorFDA)
names(res)[2]<-"zscore" # rename to zscore because some functions expect this naming
res$ae_descr="CT"
reflist_old<-res[res$drug=="ceritinib",]

# combine reflist data
reflist<-rbind(reflist_ext[,-c(3,4)],reflist_old)
reflist$drug[4]<-"ceritinib_ext"

# load raw dataa
load("rdata/data_rawdata_to_foldchange_ExtraDrugs.Rdata") # load source data to be used
alldat_fc<-alldat_fc2
alldat_fcClean <- alldat_fc[,c("gene","drug","human","logFC")] # remove pvalues, fdr etc
alldat_fcClean$drug<-tolower(alldat_fcClean$drug)
alldat_fc<-alldat_fcClean[alldat_fcClean$gene %in% signature,]

# calc mean fold change in dataset
logfc_mean_internal <- select(alldat_fc,-human) %>% 
  group_by(drug, gene) %>%
  summarise_each(funs(mean))
logfc_mean_internal$gene<-as.character(logfc_mean_internal$gene)
log_fc_drug_matrix_internal<-dcast(formula = drug~gene, value.var="logFC", data=logfc_mean_internal)
rownames(log_fc_drug_matrix_internal)<-tolower(as.character(log_fc_drug_matrix_internal$drug)) # remove drugnames column and set as rownames
log_fc_drug_matrix_internal$drug<-NULL
deltagene<-log_fc_drug_matrix_internal[c(reflist$drug),]
deltagene$PRKRIP1<-0 # not in new DGE data

deltagene_test<-deltagene
reflist_test<-reflist

matxy<-mergeObsGeneFC(reflist = reflist_test,
                      fc_matrix = deltagene_test,  ae="CT")

predy<-as.data.frame(predict(ENET1n$finalModel, newx=matxy$x,
                             s=ENET1n$bestTune$lambda))

predy$drug<-rownames(predy)
names(predy)[1]<-"predicted"
pred_new<-merge(reflist,predy,by="drug")



# Load old external predictions
load(file="results/ALLRES image_2019-08-20.Rdata")

matxy<-mergeObsGeneFC(reflist = reflist_test, fc_matrix = deltagene_test,  ae="CT")
predy<-as.data.frame(predict(ENET1n$finalModel, newx=matxy$x[,vip1$gene,drop=FALSE],s=ENET1n$bestTune$lambda))

predy$drug<-rownames(predy)
names(predy)[1]<-"predicted"
pred_old<-merge(reflist,predy,by="drug")

# Combine new and old external predictions
pred_comb<-rbind(pred_new,pred_old)


plot(pred_comb$predicted,pred_comb$zscore, pch=19,col="blue",xlim=c(0,1.7), ylim=c(0,1.7),
     xlab="Reporting odds ratio from FAERS",
     ylab="Predicting reporting odds ratio from signature",
     main="")
text(pred_comb$predicted,pred_comb$zscore+.1,labels = toupper(substr(pred_comb$drug,1,3)))
abline(0,1,lty=2)
dev.off()

