############### Preparations
rm(list=ls())

library(gplots)
library(dplyr)
library(tidyr)
library(caret)
library(reshape2)

source("00_functions.R")

############### Data aggregration

# 1: parse gene lists for enrichment 
###  # only need to be run once  # ----- # source("analysis_parsegenelists.R")

# 2: LINCS gene expression filesto fold change and gene expresion datasets
###  # only need to be run once  # ----- # source("data_rawdata_to_foldchange.r")
###  # only need to be run once  # ----- # source("data_rawdata_to_geneExpression.R")

# 3: Create analysis datasets (filtering, consistent genes) for NIB, nonNIB, and nonNIB+NIB drugs
###  # only need to be run once  # ----- # source("data_foldchange_to_matrix.R")   #  fold-change dataset matrix (for ENET)

# FAERS risk scores ; to order based on risk score drugs
load("FAERS_external/riskscoresMeta.Rdata") # generated in /data/lincs_predictCT/estimateCTrisk_finalJuly11.R
res<-resMeta
res$drug<-as.character(res$drug)
res<-select(res,drug,rorFDA)
names(res)[2]<-"zscore" # rename to zscore because some functions expect this naming
res$ae_descr="CT"
reflist<-res


########################################################################################
# Processing of data

# load raw data
load("rdata/data_rawdata_to_foldchange.Rdata") # load source data to be used
alldat_fc1<-alldat_fc #unchanged version
alldat_fcClean <- alldat_fc[,c("gene","drug","human","logFC")] # remove pvalues, fdr etc


#####



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

# calc mean fold change in dataset
logfc_mean_internal <- select(alldat_fc,-human) %>% 
  group_by(drug, gene) %>%
  summarise_each(funs(mean))
logfc_mean_internal$gene<-as.character(logfc_mean_internal$gene)
log_fc_drug_matrix_internal<-dcast(formula = drug~gene, value.var="logFC", data=logfc_mean_internal)
rownames(log_fc_drug_matrix_internal)<-tolower(as.character(log_fc_drug_matrix_internal$drug)) # remove drugnames column and set as rownames
log_fc_drug_matrix_internal$drug<-NULL
deltagene<-log_fc_drug_matrix_internal[reflist$drug,]

# calc mean pvalue in dataset
alldat_fc1 <- alldat_fc1[,c("gene","drug","human","pvalue")] # remove pvalues, fdr etc
alldat_fcP<-alldat_fc1[alldat_fc1$gene %in% genes_nib_intersect,]
logfc_mean_internal <- select(alldat_fcP,-human) %>% 
  group_by(drug, gene) %>%
  summarise_each(funs(mean))
logfc_mean_internal$gene<-as.character(logfc_mean_internal$gene)
logfc_mean_internal$drug<-tolower(logfc_mean_internal$drug)

########################################################################################
# PCA plots

pca<-prcomp(deltagene, scale=T)
tol21rainbow= c("#084594","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## Figure 1D
library(rgl)
plot3d(pca$x[,1:3], col=tol21rainbow, size=7)
text3d(pca$x[,1:3], text=toupper(substr(rownames(pca$x),start=1,stop=3)), adj=1, cex=1,col=tol21rainbow) 
rgl.snapshot("results/figure_3dpca.png",fmt="png",top=TRUE)





########################################################################################
# Jaccard index plots for top 100-500 genes 
datP<-logfc_mean_internal
top250<-data.frame()
for(d in unique(datP$drug)){
  xx<-datP[datP$drug==d,]
  xx<-as.data.frame(xx[order(xx$pvalue),][1:250,])
  top250<-rbind(xx,top250)
}

#Figure 1C - jaccard index GE data
ji250<-compJI(top250)
library(gplots)
win.metafile("results/Figure 1C.emf",width=15,height=15,pointsize = 10)
par(mar=c(5, 2, 3, 1) + 0.1, oma=c(10,5,6,7))
heatmap.2(x=ji250, 
          col =  rainbow(25,start = 0,end=.9),
          scale="none",
          margins=c(7.5,1), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, symbreaks=FALSE, 
          dendrogram='none',
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          cexRow=3, cexCol=3,
          keysize=.2, key.par=list(mar=c(3.5,0,3,0)),
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(.6, 5), lwid=c(.4, 10, 1))

dev.off()


#########################################################################
# Enrichment of raw gene signatures

enrichmentPathways<-function(genelist, explist, pval=0.01, allgenes=23896){
  allgene_in_data <- length(explist)
  res <-  lapply(1:length(genelist), FUN = enrichmentPathwaysSingle ,
                 genelist=genelist,explist=explist,allgenes=allgenes)
  res<-do.call("rbind", res)
  #res<-res[res$pval<pval,]
  res<-res[order(res$pval),]
  return(res)
}

load("rdata/data_genelists.Rdata")


###########################################################################
### Figure kinases

sum<-data.frame()
for(d in unique(top250$drug)){
  print(d)
  enr.raw<-enrichmentPathways(genelist=genelist.KEA,
                              explist=as.character(top250$gene[top250$drug==d]))
  enr.raw<-enr.raw[enr.raw$pval<.05,]
  enr.raw$drug<-d
  sum<-rbind(sum,enr.raw)
}
sumold<-sum
# remove diseases from gene lists
#sum<-sum[-grep("melanoma|measles|addiction|malaria|trypanosomiasis|rheumatoid|phagocytosis|Oocyte|leishmaniasis|platelet|influenza|anemia|sclerosis|bacterial|leukemia|diabetic|depression|amoebiasis|pertussis|legeionellosis|shigellosis|carcinogenesis|carcinoma|hepatitis|tuberculosis|glioma|toxoplasmosis|disease|cancer|infection",sum$pathway,ignore.case = T),]
#sum$ml10p<--log10(sum$pval)

res.corpw<-data.frame()
for(pw in unique(sum$pathway)){
  
  x.pval<-sum[sum$pathway==pw,c("pval","drug")]
  x.pval$dich<-1
  x.pval$mlp<--log10(x.pval$pval)
  
  mis<-setdiff(as.character( unique(sum$drug)), sum$drug[sum$pathway==pw] )
  x.pval<-rbind(x.pval, data.frame(pval=1,drug=mis,dich=0, mlp=0))
  
  xy<-merge(x.pval,reflist,by="drug")
  r<-cor(xy$mlp,xy$zscore)
  #plot(xy$mlp,xy$zscore)
  res.corpw<-rbind(res.corpw,data.frame(pw=pw,r=r))
}

res.corpw<-res.corpw[order(res.corpw$r),]
res.corpw<-res.corpw[abs(res.corpw$r)>.25,]
res.corpw<-res.corpw[order(res.corpw$r),]
res.corpw$pw<-factor(res.corpw$pw,levels=res.corpw$pw)


# Figure 3B
pKIN<-ggplot(data=res.corpw)+
  geom_bar(aes(x=pw,y=r,fill=r),stat="identity")+
  coord_flip()+
  theme_bw()+
  ggtitle("Enriched kinases")+
  scale_fill_gradientn(colours=bluered(10),guide="none")+
  ylab("Correlation (r) with cardiotoxicity risk score")+
  xlab("")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




###########################################################################
### Figure KEGG processes

sum<-data.frame()
for(d in unique(top250$drug)){
  print(d)
  enr.raw<-enrichmentPathways(genelist=genelist.KEGG,
                              explist=as.character(top250$gene[top250$drug==d]))
  enr.raw<-enr.raw[enr.raw$pval<.05,]
  enr.raw$drug<-d
  sum<-rbind(sum,enr.raw)
}
sumold<-sum
# remove diseases from gene lists
sum<-sum[-grep("melanoma|measles|alcoholism|addiction|malaria|trypanosomiasis|rheumatoid|phagocytosis|Oocyte|leishmaniasis|platelet|influenza|anemia|sclerosis|bacterial|leukemia|diabetic|depression|amoebiasis|pertussis|legeionellosis|shigellosis|carcinogenesis|carcinoma|hepatitis|tuberculosis|glioma|toxoplasmosis|disease|cancer|infection",sum$pathway,ignore.case = T),]
sum$ml10p<--log10(sum$pval)

res.corpw<-data.frame()
for(pw in unique(sum$pathway)){
  
  x.pval<-sum[sum$pathway==pw,c("pval","drug")]
  x.pval$dich<-1
  x.pval$mlp<--log10(x.pval$pval)
  
  mis<-setdiff(as.character( unique(sum$drug)), sum$drug[sum$pathway==pw] )
  x.pval<-rbind(x.pval, data.frame(pval=1,drug=mis,dich=0, mlp=0))
  
  xy<-merge(x.pval,reflist,by="drug")
  r<-cor(xy$pval,xy$zscore)
  #plot(xy$mlp,xy$zscore)
  res.corpw<-rbind(res.corpw,data.frame(pw=pw,r=r))
}
res.corpw$pw<-sapply(strsplit(as.character(res.corpw$pw),split = "_"),function(x) x[[1]])
#res.corpw<-res.corpw[-grep("mouse",res.corpw$pw),]
#res.corpw$pw<-gsub(" \\(human|\\)","",res.corpw$pw)
res.corpw<-res.corpw[order(res.corpw$r),]
res.corpw<-res.corpw[abs(res.corpw$r)>.25,]
res.corpw<-res.corpw[order(res.corpw$r),]
res.corpw$pw<-factor(res.corpw$pw,levels=res.corpw$pw)

#Figure 3D
pKEGG<-ggplot(data=res.corpw)+
  geom_bar(aes(x=pw,y=r,fill=r),stat="identity")+
  coord_flip()+
  theme_bw()+
  ggtitle("KEGG enriched terms")+
  scale_fill_gradientn(colours=bluered(10), guide = "none")+
  xlab("")+ylab("Correlation (r) with cardiotoxicity risk score")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
print(pKEGG)

library(gridExtra)

# Figure 3B 3C
pdf("results/figure_enrichment raw data for publication.pdf",height=8,width=16)
grid.arrange(pKIN, pKEGG, ncol=2)
dev.off()
























