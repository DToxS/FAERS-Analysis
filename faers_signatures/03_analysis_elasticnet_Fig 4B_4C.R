rm(list=ls())
library(caret)
library(doParallel)
library(parallel)
library(dplyr)
library(reshape2)

source("functions.R")

# FAERS risk scores 
load("FAERS_external/riskscoresMetaOrig.Rdata")
res<-resMeta
res$drug<-as.character(res$drug)
res<-select(res,drug,rorFDA)
names(res)[2]<-"zscore" # rename to zscore because some functions expect this naming
res$ae_descr="CT"
reflist<-res

# load raw dataa
load("rdata/data_rawdata_to_foldchange.Rdata") # load source data to be used

alldat_fcClean <- alldat_fc[,c("gene","drug","human","logFC")] # remove pvalues, fdr etc

#  alldat_fcClean$drug<-tolower(alldat_fcClean$drug)
#  alldat_fcClean<-alldat_fcClean[alldat_fcClean$drug %in% reflist$drug,]
  

# filter for genes present in all cell lines
genes_nib_intersect<-character(0); c=1
for(celline in c("A","B","D","E")){
  print(celline)
  for(drug in unique(alldat_fcClean$drug[alldat_fcClean$human==celline])){
    print(length(genes_nib_intersect))
    genes<-unique(c(as.character(alldat_fcClean$gene[alldat_fcClean$drug == drug & alldat_fcClean$human==celline])))
    if(c==1) {genes_nib_intersect <- genes } else {
      genes_nib_intersect<-intersect(genes_nib_intersect,
                                     genes)
    }
    c=c+1
  }
}
alldat_fc<-alldat_fcClean[alldat_fcClean$gene %in% genes_nib_intersect,]

#alldat_fc<-alldat_fcClean[alldat_fcClean$gene %in% vip1$gene,]

# calc mean fold change in dataset
logfc_mean_internal <- select(alldat_fc,-human) %>% 
  group_by(drug, gene) %>%
  summarise_each(funs(mean))
logfc_mean_internal$gene<-as.character(logfc_mean_internal$gene)
log_fc_drug_matrix_internal<-dcast(formula = drug~gene, value.var="logFC", data=logfc_mean_internal)
rownames(log_fc_drug_matrix_internal)<-tolower(as.character(log_fc_drug_matrix_internal$drug)) # remove drugnames column and set as rownames
log_fc_drug_matrix_internal$drug<-NULL
deltagene<-log_fc_drug_matrix_internal[reflist$drug,]

# select test drugs
extdrug=reflist$drug[c(4,18)]

# settings bootstrap and cross validation
nbs=1000
nfold=5

# start a cluster
stopCluster(cl)
#nodelist <- sort( rep(c(paste("node0",1:9,sep=""),paste("node",10:16,sep="")),7))
nodelist<-rep("localhost",62)
cl<-makePSOCKcluster(nodelist)
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
clusterEvalQ(cl, library(caret))
clusterEvalQ(cl, library(doParallel))
clusterEvalQ(cl, setwd("/home/coen/faerslincs/"))


# Setup gene Expression datadata for regression
deltagene_train<-deltagene[-which(rownames(deltagene) %in% extdrug),]
reflist_train<-reflist[-which(reflist$drug %in% extdrug),]
deltagene_test<-deltagene[which(rownames(deltagene) %in% extdrug),]
reflist_test<-reflist[which(reflist$drug %in% extdrug),]

matxy<-mergeObsGeneFC(reflist = reflist_train, fc_matrix = deltagene_train, ae="CT")
x1=matxy$x
y1=matxy$y


# create seeds 
set.seed(123)
seed1<-lapply(1:1001,function(x)sample(100:10000,size=60,replace=FALSE))
egrid <- expand.grid(.alpha =exp(seq(-6,0,length=40)),.lambda = exp(seq(-4,4,length=40)))

set.seed(123)
bsdat<-t(sapply(1:nbs, function(x) sample(1:length(y1), size = length(y1), replace=T)))


# run bootstrapped elastic net
allres<-data.frame()
allresbs<-data.frame()
count=1

rr<-data.frame()
for(i in 1:nbs){
  cat(" ")
  cat(i)
  set.seed(12345)
  ENET1 <- train(x=as.matrix(x1[bsdat[i,],]), 
                 y=y1[bsdat[i,]],
                 method = "glmnet",
                 tuneGrid = egrid,
                 trControl=trainControl(method='repeatedcv',number=nfold,repeats=30,seeds=seed1))
  
  vip1<-vip_model(ENET1)
  vip1<-  vip1[  vip1$Overall!=0 & !is.nan(vip1$Overall),]
  
  if(nrow(vip1)>0) vip1$id=i
  if(nrow(vip1)>0) vip1$ae='CT'
  if(nrow(vip1)>0) rr<-rbind(rr, vip1)
}
cat("\n")
print("finish bs ")
# do final fit
sg<-sort(  table(rr$gene), decreasing = T)
sg<-as.data.frame(sg)
names(sg)[1]<-"gene"
meanval<-tapply(rr$Overall, rr$gene, mean)
mval<-data.frame(names(meanval), meanval)
names(mval)[1]<-"gene"
mval$gene<-as.character(mval$gene)
mval<-mval[order(mval$gene),]
sg$gene<-as.character(sg$gene)
sg<-sg[order(sg$gene),]
sg2<-cbind(sg, mval)
sg2$Freq<-(sg2$Freq/nbs)*100
sg2$freqmeanrat<-sg2$meanval*sg2$Freq
allres<-sg2
allresbs<-rr

#  selected genes in bootstrap
res_ct<-allres
res_ct<-res_ct[order(-res_ct$Freq),]
print("run enet for quantiles")
#  run final elastic net accross different percentiles
r2rmse<-data.frame()
ae_i="CT"
for(q in c(seq(.95,.99,.01),seq(.99,0.9995,.00005))){
  #for(q in c(seq(.98005,0.9995,.0005))){
  cat(" ") ; cat(q)
  matxy<-mergeObsGeneFC(reflist = reflist_train,
                        fc_matrix = deltagene_train, ae="CT")
  x=matxy$x
  y=matxy$y
  
  # edit june
  
  #select<-res_ct[res_ct$Freq>=quantile(res_ct$Freq, q,na.rm=T) &                  !is.na(res_ct$Freq),]
  select<-res_ct[res_ct$freqmeanrat>=quantile(res_ct$freqmeanrat, q,na.rm=T) &                  !is.na(res_ct$freqmeanrat),]
  if(nrow(select)>1){
    set.seed(123)
    ENET1  <- train(x=as.matrix(x1[,select$gene]), y=y1, 
                    method = "glmnet",
                    tuneGrid = egrid,
                    trControl=trainControl(method='repeatedcv',number=nfold,repeats=30,seeds=seed1))
    
    res<-ENET1$results[ENET1$results$alpha==ENET1$best$alpha &ENET1$results$lambda==ENET1$best$lambda,]
    res$ae=ae_i
    res$q=q 
    
    
    #vip1<-res_ct[res_ct$Freq>=quantile(res_ct$Freq, q,na.rm=T),]
    vip1<-res_ct[res_ct$freqmeanrat>=quantile(res_ct$freqmeanrat, q,na.rm=T),]
    
    ENET1n  <- train(x=as.matrix(x1[,vip1$gene[!is.na(vip1$gene)]]), y=y1,
                     tuneGrid = egrid, method = "glmnet",
                     trControl=trainControl(method='repeatedcv',number=nfold,repeats=50,seeds=seed1))
    
    matxy<-mergeObsGeneFC(reflist = reflist_test, fc_matrix = deltagene_test,  ae="CT")
    predy<-predict(ENET1n$finalModel, newx=matxy$x[,vip1$gene,drop=FALSE],s=ENET1n$bestTune$lambda)
    
    res$reg=predy[1,1]
    res$sun=predy[2,1]
    
    r2rmse<-rbind(r2rmse,res)
  }
  
  
}






minRMSE<-lapply(unique(r2rmse$ae), function(a) {
  aa<-r2rmse[r2rmse$ae==a,]
  return(aa[aa$RMSE==min(aa$RMSE),])
})
minRMSE<-do.call("rbind",minRMSE)

minR2<-lapply(unique(r2rmse$ae), function(a) {
  aa<-r2rmse[r2rmse$ae==a,]
  return(aa[aa$Rsquared==max(aa$Rsquared),])
})
minR2<-do.call("rbind",minR2)


# select based on optimal percentile
ae_i="CT"
q=minRMSE$q[minRMSE$ae==ae_i][nrow(minRMSE)]



vip1<-res_ct[res_ct$freqmeanrat>=quantile(res_ct$freqmeanrat, q,na.rm=T),]
ENET1n  <- train(x=as.matrix(x1[,vip1$gene[!is.na(vip1$gene)]]), y=y1,
                 tuneGrid = egrid, method = "glmnet",
                 trControl=trainControl(method='repeatedcv',savePredictions = TRUE,number=nfold,repeats=50,seeds=seed1))
predy<-predict(ENET1n$finalModel, newx=matxy$x[,vip1$gene,drop=FALSE],s=ENET1n$bestTune$lambda)

matxy<-mergeObsGeneFC(reflist = reflist_test, fc_matrix = deltagene_test,  ae="CT")
predy<-predict(ENET1n$finalModel, newx=matxy$x[,vip1$gene,drop=FALSE],s=ENET1n$bestTune$lambda)
genes.selected<-vip_model(ENET1n)

save(r2rmse,allresbs,allres,vip1,predy,ENET1n,genes.selected,file=paste("results/ALLRES_2019-08-20_",paste(extdrug,collapse="-"),"nfold",nfold,"nbs",nbs,"ndesc",nrow(vip1),".Rdata",sep="_"))
save.image(file="results/ALLRES image_2019-08-20.Rdata")

load(file="results/ALLRES_2019-08-20__regorafenib-sunitinib_nfold_5_nbs_1000_ndesc_26_.Rdata")
load(file="results/ALLRES image_2019-08-20.Rdata")


## Figure external prediction observed 

matxy<-mergeObsGeneFC(reflist = reflist_test, fc_matrix = deltagene_test,  ae="CT")
predy<-predict(ENET1n$finalModel, newx=matxy$x[,vip1$gene,drop=FALSE],s=ENET1n$bestTune$lambda)




pdf(file="figextpredold.pdf",height=6,width=6)
plot(matxy$y,predy, pch=19,col="blue",xlim=c(0,2), ylim=c(0,2),
     xlab="Reporting odds ratio from FAERS",
     ylab="Predicting reporting odds ratio from signature",
     main="")
text(matxy$y,predy+.1,labels = c("REG","SUN"))
abline(0,1,lty=2)
 dev.off()

# Figure importance final signature
library(caret)
library(glmnet)
vipFinal<-vip_model(ENET1n)
vipFinal<-vipFinal[order(vipFinal$Overall),]
vipFinal<-vipFinal[vipFinal$Overall!=0,]
vipFinal$gene<-factor(vipFinal$gene,levels=vipFinal$gene)


############## Figure 4B

pdf("results/Figure 4B.pdf",width=6,height=6)
ggplot(data=vipFinal)+
  geom_bar(aes(y=Overall,x=gene,fill=Overall),stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_gradient(guide="none")+
  xlab("Genes")+ylab("Relative importance (%)")
dev.off()




# Figure importance of each drug
ENET1<-ENET1n
out_of_fold <- ENET1$pred
out_of_fold <- out_of_fold[out_of_fold$alpha==ENET1$bestTune$alpha & out_of_fold$lambda==ENET1$bestTune$lambda,]
out_of_fold$drug<-rownames(ENET1$training)[out_of_fold$rowIndex]
out_of_fold$re<-abs(out_of_fold$pred-out_of_fold$obs)/(abs(out_of_fold$pred)+abs(out_of_fold$obs))
inf<-aggregate(out_of_fold$re, by=list(out_of_fold$drug), mean)
inf$pred<-aggregate(out_of_fold$pred, by=list(out_of_fold$drug), mean)[,2]
inf$obs<-aggregate(out_of_fold$obs, by=list(out_of_fold$drug), mean)[,2]
inf<-inf[order(inf$x),]
names(inf)<-c("drug","meanAbsPE","pred","obs")
inf$RE<-inf$pred-inf$obs
inf$meanAbsPE<-inf$meanAbsPE*100
inf<-inf[order(inf$meanAbsPE),]
inf$drug<-factor(inf$drug,levels=inf$drug)

pdf("results/figure_importance of drugs for publication.pdf",width=6,height=6)
ggplot(data=inf,aes(y=meanAbsPE,x=drug))+
  geom_bar(stat="identity")+coord_flip()+
  ylab("Symmetric mean absolute percentage error (%)")+xlab("Kinase inhibitor")+
  theme_bw()
dev.off()

# Figure obs-pred training


sum<-cbind(aggregate(out_of_fold$obs, by=list(drug=out_of_fold$drug), mean), 
           aggregate(out_of_fold$pred, by=list(drug=out_of_fold$drug), mean)[,2], 
           aggregate(out_of_fold$pred, by=list(drug=out_of_fold$drug), sd)[,2])
names(sum)<-c("drug","obs","mean","sd")
sum$ae<-"CT"

ct<-ENET1$results[ENET1$results$alpha==ENET1$best$alpha &ENET1$results$lambda==ENET1$best$lambda,]
ct$ae<-"CT"
resbest<-ct
resbest$x=0.4
resbest$y=1.5
resbest$labelR2=paste("R^2==",round(resbest$Rsquared,2))
resbest$labelRMSE=paste("RMSE==",round(resbest$RMSE,3))


sum$err<-sum$obs-sum$mean
sum$label<-sum$drug
sum$label[abs(sum$err)<.1]<-NA

##### Figure 4C

pdf(file='Figure 4C.pdf',width=6, height=6)
ggplot(data=sum)+
  geom_point(aes(x=obs,y=mean,group=drug,col=drug))+
  geom_text(aes(x=obs,y=mean+.05,label=label,group=drug,col=drug))+
  geom_errorbar(aes(x=obs,ymin=mean-sd,ymax=mean+sd,group=drug,col=drug))+
  geom_abline(intercept=0)+theme_bw()+xlab("Observed risk score")+ylab("Predicted risk score")+
  geom_text(data=resbest,aes(x=x,y=y,label=labelRMSE), parse=T)+
  geom_text(data=resbest,aes(x=x,y=y-.1,label=labelR2), parse=T)+
  theme(legend.position="bottom")+
  coord_cartesian(xlim=c(-0,2),ylim=c(0,2))+
  scale_color_discrete(guide="none")

dev.off()

#### Figure heatmap signature DEG

x_ct <- mergeObsGeneFC(reflist = reflist, fc_matrix = deltagene, ae="CT")$x

library(gplots)
library(gridGraphics)
library(gridExtra)
vip1<-vip_model(ENET1n)
vip1<-vip1[vip1$Overall>0,]
vip1<-vip1[order(vip1$Overall),]
x_ct<-x_ct[reflist$drug,]
x_ct<-x_ct[,vip1$gene]

pdf("results/figure Heatmap Signature.pdf")
heatmap.2(x_ct[,as.character(vip1$gene)],
                              dendrogram='none', Rowv=FALSE, Colv=FALSE,
                              col=rainbow(40),cexRow=1.2, cexCol=1.2,
                              margins=c(15,10),key=TRUE,
                              trace="none", main="",
                              key.title = "")
  dev.off()
  
grid.newpage()
pdf("figures/figure_heatmap_signature_VDHF.pdf", width=15, height=15)
grid.arrange(grobs=gl, clip=TRUE,ncol=2)
dev.off()




