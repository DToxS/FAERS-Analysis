rm(list=ls())
library(tidyr)
library(reshape2)
library(ggplot2)
library(meta)
library(xlsx)
setwd("")

tkiNames<-tolower(c("axitinib","bosutinib","cabozantinib","ceritinib","crizotinib","dabrafenib","dasatinib",
                    "erlotinib","gefitinib","imatinib","lapatinib","nilotinib","pazopanib","ponatinib",
                    "regorafenib","ruxolitinib","sorafenib","sunitinib","tofacitinib","trametinib",
                    "vandetanib","vemurafenib","afatinib"))

####################################################################################
# Selected ADRs from MEDDRA ontology

meddra<-read.csv("./reference_files/meddra_annotate.csv",stringsAsFactors = F)
med<-meddra$level4[meddra$level2 %in% c("heart failures","cardiac disorder signs and symptoms","myocardial disorders")]

####################################################################################
# FDA AERSMINE analysis
####################################################################################
dat<-read.csv("./reference_files/aersMineExploreDataSet_1957.tsv",sep="\t",skip=7,stringsAsFactors = F)
dat<-dat[,c(grep("Total_Adverse.Events_Reports|Adverse.Events",names(dat)),grep("Absolute.Counts",names(dat)))]
names(dat)<-gsub("Absolute.Counts_|___aggregated","",names(dat))
names(dat)[1:2]<-c("ADR","TotalReports")
dat<-dat[,-grep("ib.1",names(dat))]

dat2<-melt(data = dat,id.vars = c("ADR","TotalReports"))
dat2$value<-gsub(",","",dat2$value)
dat2$TotalReports<-gsub(",","",dat2$TotalReports)
dat2$value<-as.numeric(as.character(dat2$value))
dat2$TotalReports<-as.numeric(as.character(dat2$TotalReports))
dat2<-dat2[dat2$ADR!="Unique Patients",]
dat2$variable<-as.character(dat2$variable)
dat2$value[is.na(dat2$value)]<-0

# Calculate ORs
res<-data.frame()
for(d in unique(dat2$variable)){
  print(d)
  dycy<-sum(dat2$value[dat2$variable ==d & dat2$ADR %in% med])
  dycn<-sum(dat2$value[dat2$variable ==d & !dat2$ADR %in% med])
  
  dncy<-sum(dat2$TotalReports[dat2$variable !=d & dat2$ADR %in% med])
  dncn<-sum(dat2$TotalReports[dat2$variable !=d & !dat2$ADR %in% med])
  
  dncy<-dncy - dycy
  dncn<-dncn - dycn
  
  orr=(dycy / dycn)/ (dncy / dncn)
  ln_orr<-log(orr)
  SE_ln_orr <- sqrt((1/dycy) + (1/dycn) + (1/dncy) + (1/dncn))
  lnCIup<-ln_orr+1.96 * SE_ln_orr
  lnCIdown<-ln_orr-1.96 * SE_ln_orr
  CIup<-exp(lnCIup)
  CIdown<-exp(lnCIdown)
  
  res<-rbind(res,data.frame(drug=d,orr=orr,SElnORR=SE_ln_orr,CIup,CIdown,dycy=dycy,dycn=dycn))
  
}
res<-res[order(res$orr),]
res$drug<-factor(res$drug,levels=res$drug)
resFAERS<-res

####################################################################################
# WHO dataset extracted from http://www.vigiaccess.org/ on July 11 2017

# Import data from XLS sheet
library(xlsx)

vigiData<-data.frame()
for( tki in tkiNames){
  dat<-xlsx::read.xlsx("./reference_files/vigiAccessTKI.xlsx",tki,header=F)
  dat$X1<-as.character(dat$X1)
  dat$tki<-tki
  dat$ae<-sapply(strsplit(dat$X1," \\("),function(x) x[[1]])
  dat$nr<-as.numeric(gsub("\\)","",sapply(strsplit(dat$X1," \\("),function(x) x[[2]])))
  dat$X1<-NULL
  vigiData<-rbind(vigiData,dat)
}
vigiData$ae<-as.character(tolower(vigiData$ae))

# Calculate ORs
res<-data.frame()
for( tki in tkiNames){
  
  NoTki_other<-sum(vigiData$nr[ vigiData$ae =="total"])
  Tki_other<-vigiData$nr[vigiData$tki==tki & vigiData$ae =="total"]
  
  NoTki_card<-sum(vigiData$nr[vigiData$tki!=tki & vigiData$ae !="total" & vigiData$ae %in% med])
  Tki_card<-sum(vigiData$nr[vigiData$tki==tki & vigiData$ae !="Total" &  vigiData$ae %in% med])
  
  NoTki_other<-NoTki_other- NoTki_card
  Tki_other<-Tki_other -Tki_card
  
  dycy=Tki_card
  dycn=Tki_other
  dncy=NoTki_card
  dncn=NoTki_other
  orr=(dycy / dycn)/ (dncy / dncn)
  
  ln_orr<-log(orr)
  SE_ln_orr <- sqrt((1/dycy) + (1/dycn) + (1/dncy) + (1/dncn))
  lnCIup<-ln_orr+1.96 * SE_ln_orr
  lnCIdown<-ln_orr-1.96 * SE_ln_orr
  CIup<-exp(lnCIup)
  CIdown<-exp(lnCIdown)
  
  res<-rbind(res,data.frame(drug=tki,orr=orr,SElnORR=SE_ln_orr,CIup,CIdown,dycy=dycy,dycn=dycn))
  
}
res<-res[order(res$orr),]
resWHO<-res

####################################################################################
# Compare datasets

################################## Figure 1C

# Compare ROR
res2<-merge(resFAERS,resWHO,by="drug")
pdf("Figure 1C.pdf",width=9,height=9)
plot(res2$orr.x,res2$orr.y,xlim=c(0.2,3),ylim=c(0.2,3),type="n",log="xy",
     xlab="ROR FAERS",ylab="ROR WHO")
abline(0,1,lty=2)
text(res2$orr.x,res2$orr.y,labels=res2$drug)
dev.off()


##################################################################
# Meta analysis of ROR
names(res2)
res3<-res2[,c(1,2,3,4,5,8,9,10,11)]
names(res3)<-c("drug","orrFDA","seFDA","CIdFDA","CIuFDA","orrVIG","seVIG","CIdVIG","CIuVIG")

resMeta<-data.frame()
for(tki in tkiNames){

  rorFDA=res3$orrFDA[res3$drug==tki]
  rorVIG=res3$orrVIG[res3$drug==tki]
  seFDA=res3$seFDA[res3$drug==tki]
  seVIG=res3$seVIG[res3$drug==tki]
  
  or <- c(rorFDA,rorVIG)
  se <- c(seFDA,seVIG)
  logor <- log(or)
  
  or.fem <- metagen(logor, se, sm = "OR")
  metaOR<-exp(or.fem$TE.fixed)
  metaCIdown<-exp(or.fem$lower.fixed)
  metaCIup<-exp(or.fem$upper.fixed)
  resMeta<-rbind(resMeta,data.frame(drug=tki,
                                    rorFDA=rorFDA,
                                    CIdFDA=res3$CIdFDA[res3$drug==tki],
                                    CIuFDA=res3$CIuFDA[res3$drug==tki],
                                    rorVIG=rorVIG,
                                    CIdVIG=res3$CIdVIG[res3$drug==tki],
                                    CIuVIG=res3$CIuVIG[res3$drug==tki],
                                    metaOR=metaOR,
                                    metaSE=or.fem$seTE.fixed,
                                    metaCIup=metaCIup,
                                    metaCIdown=metaCIdown))
}
resMeta<-resMeta[order(resMeta$rorFDA),]
resMeta$drug<-factor(resMeta$drug,levels=resMeta$drug)

save(x=resMeta,file="./riskscoresMetaOrig.Rdata")



######################## Figure 1B

# Figure 1B
pdf("figure 1B.pdf",width=7, height=7)
ggplot(data=resMeta)+
  geom_segment(aes(x=drug,xend=drug,y=CIdFDA,yend=CIuFDA),size=1.1,col='darkgray')+
  geom_point(aes(x=drug,y=rorFDA,col=rorFDA),size=1.8)+
  coord_flip()+
  theme_bw()+
  ylab("Reporting odds ratio")+xlab("Kinase inhibitor")+
  scale_color_gradient(low="blue",high="red",guide = "none")
dev.off()


