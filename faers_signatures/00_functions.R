############################################################## ###############################  
# data_processing_GE_FC_beforeAnalysis
############################################################## ###############################  

# Determine top ranking genes for each drug
topRank<-function(keep=100, deltagene=deltagene_unfiltered){
  topgenes<-data.frame()
  for(d in unique(rownames(deltagene))){
    print(d)
    genesD<-deltagene[which(rownames(deltagene)==d),]
    top<-genesD[order(-abs(genesD))][1:keep]
    topgenes<-rbind(topgenes,
                    data.frame(d=d,gene=names(top), logFC= as.numeric(top)))
  }
  
  dg<-data.frame()
  for(g in as.character(unique(topgenes$gene))){
    di<-as.character(topgenes$d[topgenes$gene==g])
    dg<-rbind(dg,data.frame(drug=di,gene=g))
  }
  return(dg)
}

# Determine gene frequencies in dataset
countGenes<-function(data, countThreshold=10){
  counts<-data.frame(count=as.numeric(table(data$gene)), gene=names(table(data$gene)))
  counts<-counts[counts$count>=countThreshold,]
  counts<-counts[order(-counts$count),]
  data<-data[data$gene %in% counts$gene,]
  data$gene<-factor(as.character(data$gene), levels=unique(counts$gene))
  data$drug<-data$d
  data$gene<-data$gene
  return(data)
}

#Convert factor to integer
factor.to.int <- function(f) {
  (as.integer(f) - 1) / (length(levels(f)) - 1)
}

#Calculate JI from top ranking genes
compJI<-function(dg250){
  jiall<-matrix(nrow=length(unique(dg250$drug)),ncol=length(unique(dg250$drug)))
  c=1
  for(d in unique(dg250$drug)){
    genesi<-as.character(dg250$gene[dg250$drug==d])
    
    ji<-sapply(unique(dg250$drug), function(i) {
      genes_drug_i<-as.character(dg250$gene[dg250$drug==i])
      inters<-length(intersect(genes_drug_i,genesi))
      unions<-length(unique(c(genes_drug_i,genesi)))
      ji<-inters/unions
      return(ji)
    })
    ji<-as.matrix(ji)
    rownames(ji)<-unique(dg250$drug)
    colnames(ji)<-d
    jiall[,c]<-ji
    c=c+1
  }
  rownames(jiall)<-unique(dg250$drug)
  colnames(jiall)<-unique(dg250$drug)
  return(jiall)
}


# function used in the final regression analysis aug 2017
mergeObsGeneFC2<-function(reflist, fc_matrix, ae){
  # select AE 
  
  reflist_i <- reflist[reflist$ae_descr == ae,]
  
  fc_matrix$drug<-rownames(fc_matrix)
  fc_matrix<-fc_matrix[,c(ncol(fc_matrix),1:(ncol(fc_matrix)-1))]
  
  # merge datasets to matrix
  mat<-merge(reflist_i, fc_matrix, by="drug")
  rownames(mat)<-mat$drug
  mat$drug<-NULL
  mat$ae_descr<-NULL
  y<-mat$zscore;   mat$zscore<-NULL
  mat<-as.matrix(mat)
  isnan<-which(is.nan(y))
  if(length(isnan)>0){
    mat<-mat[-isnan,]
    y<-y[-isnan]
  } # remove nan records if needed
  # return descriptor matrix and z-scores
  return(list(x=mat,y=y))
}


# Function to match genes in GE and FC datasets
matchGenesGE_FC<-function(dataGE, dataFC){
  genesGE<-colnames(dataGE)
  genesFC<-colnames(dataFC)
  genes_select<-intersect(genesFC,genesGE)
  dataGE2<-dataGE[,genes_select]
  dataFC2<-dataFC[,genes_select]
  return(list(dataGE=dataGE2, dataFC=dataFC2))
}

# Function to remove lowly expressed genes
filterData<-function(GE_FC, GEfilterPerc, GEfilterThres){ 
  #GEfilterThres= count threshold
  #GEfilterPerc = percentage of samples belowThreshold
  if(GEfilterThres[1]==0 & GEfilterPerc[1] ==0 ){ 
  } else {
    belowThres<-apply(GE_FC$dataGE,2, function(x) round(length(x[x<GEfilterThres])/length(x),4)) 
    remove_genes<-colnames(GE_FC$dataGE)[which(belowThres>GEfilterPerc)]  
    if(length(remove_genes)>0){
      remove_cols<-which(colnames(GE_FC$dataGE) %in% remove_genes)  
      GE_FC$dataGE<-GE_FC$dataGE[,-remove_cols] 
      GE_FC$dataFC<-GE_FC$dataFC[,-remove_cols] 
    }
  } 
  return(GE_FC) 
}

#2-log transform + 1 all GE data
logGE<-function(GE_FC,GElog){
  if(GElog==1){
    GE_FC$dataGE<-log2(GE_FC$dataGE+1)
  }
  return(GE_FC)
}

########################################################################################################################
# analysis WGCNA data matrix to clusters
########################################################################################################################


# FUNCTION to identify optimal power
findPower.WGCNA<-function(powers = c(3:8), data){
  sft       <- pickSoftThreshold(data, powerVector = powers, verbose = 5)
  softPower <- sft$powerEstimate+1
  x         <- goodSamplesGenes(data, verbose = 3) # check if samples are OK
  if(!x$allOK) warning("check data")
  return(softPower)
}


# Function to generate the clusters
runClust.WGCNA<-function(minModSize, data, pwr, colnames, deepSplit){
  
  adjacency     <- adjacency(data, power=pwr, type="signed") # adjency matrix
  TOM           <- TOMsimilarity(adjacency, TOMType="signed") #topological overlap matrix 
  dissTOM       <- 1-TOM
  geneTree      <- fastcluster::hclust(as.dist(dissTOM), method="average")# clustered gene tree
  dynamicMods   <- cutreeDynamic(dendro            = geneTree, #generating modules 
                                 distM             = dissTOM, 
                                 deepSplit         = deepSplit,
                                 pamRespectsDendro = FALSE,
                                 minClusterSize    = minModSize)
  dynamicColors <- labels2colors(dynamicMods)
  
  # Calculate eigengenes
  MEList = moduleEigengenes(data, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs); # Calculate dissimilarity of module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average"); # Cluster module eigengenes
  
  res           <- data.frame(modules  = dynamicColors,
                              genes    = colnames(data))
  
  return(list(adjacency   = adjacency,
              TOM   = TOM, # used for visualization of networks
              geneTree    = geneTree,
              pwr=pwr,
              dynamicMods = dynamicMods,
              dynamicColors = dynamicColors,
              ME = list(colnames=colnames, ME=MEs, MEDiss=MEDiss, METree=METree),
              results     = res))
  
}

########################################################################
# WGCNA clusters to metrics
########################################################################


# Function: calculate sumMetrics for each cluster - for controls
calcSumMetrics.control<-function(control_i, cluster.nib){
  sel     <-cluster.nib$results
  
  # summary metrics CONTROLS for cluster M  
  sumMetrics_control<-data.frame()
  for(ss in unique(sel$modules)){
    print(ss)
    genes_i    <-sel$genes[sel$modules %in% ss]
    selectcols<-which(colnames(control_i) %in% genes_i)
    control_i_c<-control_i[,selectcols]
    mean_i     <-apply(control_i_c,2,mean)
    med_i      <-apply(control_i_c,2,median)
    sumMetrics_control<-rbind(sumMetrics_control, data.frame(module=ss, drug="CONTROL", 
                                                             nclus=ncol(control_i_c), mean_i, med_i))
  }
  return(sumMetrics_control)
}

# Function: calculate sumMetrics for each cluster - for drugs
# this function gets called by parLapply each time feeding a new cluster (= ss)
calcByModuleDrugsMetrics<-function(ss,sel, drugs, drugs_i){
  genes_i    <-sel$genes[sel$modules %in% ss] # sel = gene expression dataset
  drugs_i_c  <-drugs_i[,which(colnames(drugs_i) %in% c("human","plate","drug",as.character(genes_i)))]
  
  sumMetrics_drugs_i<-data.frame()
  for(dd in drugs){    # summary metrics accross replicates
    drugs_i_c_d<-    drugs_i_c[drugs_i_c$drug==dd, -c(1:3)]
    mean_i     <-apply(drugs_i_c_d,2,mean)
    #med_i      <-apply(drugs_i_c_d,2,median)
    sumMetrics_drugs_i<-rbind(sumMetrics_drugs_i, 
                              data.frame(module=ss, drug=dd, nclus=ncol(drugs_i_c_d),
                                         mean_i))
  }
  return(sumMetrics_drugs_i)
}


# combine eigengenes and drug names for calcDeltaMetrics()
calcMEvals<-function(cluster){
  library(reshape2)
  MEdrug<-cbind( cluster$ME$colnames, cluster$ME$ME )
  names(MEdrug)[1]<-"drug"
  MEdrug2<-melt(data=MEdrug, id.vars="drug")
  names(MEdrug2)<-c("drug","ME","value")
  out<-aggregate(x=MEdrug2$value,by=list(drug=MEdrug2$drug,ME=MEdrug2$ME), mean)
  out$ME<-gsub(pattern = "ME",replacement = "", out$ME)
  return(out)
}


### Function: calculate difference between controls and drug values
calcDeltaMetrics<-function(ss, ME, control_d, select_d){
  
  # start generating the predictor dataset
  deltaMetrics<-data.frame()
  control_i<-control_d[control_d$module==ss,]
  
  ME_control_i<-ME[ME$ME == ss & ME$drug== "CONTROL",]
  for(dd in unique(select_d$drug)){
    select_i<-select_d[select_d$drug==dd & select_d$module == ss, ]
    ME_i<-ME[ME$ME == ss & ME$drug== dd,]
    
    # fold change values based on total gene expression metric drug vs control
    deltaLFCsum    <- log2(sum(select_i$mean_i)/      sum(control_i$mean_i))
    #deltaLFCmax    <- log2(max(select_i$mean_i)/      max(control_i$mean_i))
    #deltaLFCmedian <- log2(median(select_i$mean_i)/   median(control_i$mean_i))
    deltaLFCmean   <- log2(mean(select_i$mean_i)/     mean(control_i$mean_i))
    
    # eigen genes
    deltaME        <- ME_i$x - ME_control_i$x 
    deltaMEnd      <- ME_i$x 
    
    deltaMetrics<-rbind(deltaMetrics, 
                        data.frame(module=ss,
                                   drug=dd,
                                   nmod=nrow(select_i),
                                   deltaLFCsum   =deltaLFCsum   ,
                                   #  deltaLFCmax   =deltaLFCmax   ,
                                   #   deltaLFCmedian=deltaLFCmedian,
                                   deltaLFCmean  =deltaLFCmean,
                                   deltaME      =deltaME      ,
                                   deltaMEnd    =deltaMEnd    ))
  }
  return(deltaMetrics)
}


# Function: calculate sumMetrics by cluster for each drug (using function above)
calcSumMetrics.drugs<-function(drugs_i, cluster.nib,drugs,cl=cl){
  sel     <-cluster.nib$results
  
  # calculate sum metrics for drugs
  sumMetrics_drugs <- data.frame()
  out<-  parLapply(cl=cl, 
                   X=unique(sel$modules),
                   fun=calcByModuleDrugsMetrics,
                   sel=sel, drugs=drugs, drugs_i=drugs_i)
  sumMetrics_drugs<-do.call("rbind",out)
  return(sumMetrics_drugs)
}



####################################################################################################
# analysis_regression.R
####################################################################################################


# Function to clean the FAERS Z-score data file prior to processing
cleanReflist<-function(reflist){
  reflist<-reflist[!is.nan(reflist$zscore),]                 # only keep available z-scores
  remove<-c("CAD","MI","ARRHYTH","QT_TDP","MYAC_NECROSIS")   # only keep relevant CT types
  reflist<-reflist[!reflist$ae_descr %in% remove,]           
  reflist<-reflist[grep("nib|trastuzumab",reflist$drug), ]# remove non-TKIs and TKIs for which no data exists
  return(reflist)
}

# Function to merge FAERS Z-score data with cluster-metrics data for regression
mergeObsClust<-function(reflist, deltaMetrics, ae, metric){
  
  # select AE and metric
  reflist_i <- reflist[reflist$ae_descr == ae,]
  deltaMetrics_m <- deltaMetrics[,c(1:3,which(names(deltaMetrics) == metric))]
  names(deltaMetrics_m)[grep("delta",colnames(deltaMetrics_m))]<-"delta"
  
  # merge datasets to matrix
  dat_i<-merge(reflist_i, deltaMetrics_m, by="drug")
  mat<-reshape2::dcast(dat_i, zscore+drug~module, value.var="delta")
  rownames(mat)<-mat$drug;   mat$drug<-NULL
  y<-mat$zscore;   mat$zscore<-NULL
  mat<-as.matrix(mat)
  return(list(x=mat,y=y))
}

# Function to create X-Y fold change matrix and matching observations
mergeObsGeneFC<-function(reflist, fc_matrix, ae){
  # select AE 
  reflist_i <- reflist[reflist$ae_descr == ae,]
  reflist_i$orr<-NULL
  
  fc_matrix$drug<-rownames(fc_matrix)
  fc_matrix<-fc_matrix[,c(ncol(fc_matrix),1:(ncol(fc_matrix)-1))]
  
  # merge datasets to matrix
  mat<-merge(reflist_i, fc_matrix, by="drug")
  rownames(mat)<-mat$drug
  mat$drug<-NULL
  mat$ae_descr<-NULL
  y<-mat$zscore;   mat$zscore<-NULL
  mat<-as.matrix(mat)
  isnan<-which(is.nan(y))
  if(length(isnan)>0){
    mat<-mat[-isnan,]
    y<-y[-isnan]
    } # remove nan records if needed
  # return descriptor matrix and z-scores
  return(list(x=mat,y=y))
}

# Function to obtain the variance importance from a regression object
vip_model<-function(fit){
  library(caret)
  vip_model<-varImp(fit)
  vip_model<-vip_model$importance
  vip_model$gene<-rownames(vip_model)
  vip_model[order(-vip_model$Overall),]
  return(vip_model)
}


# Function to perform ENET regression on cluster objects.
fitREG_LOO<-function(x,y,matxy){
  #  x=matxy$x; y=matxy$y
  
  
  egrid<-expand.grid(.alpha = seq(.001, 1, length = 40),  .lambda = seq(0.01, 1, length=40))
  ENET1  <- train(x=as.matrix(x), y=y,tuneGrid = egrid, 
                  #method = "glmnet",trControl=trainControl(method='LOOCV'))
                  method = "glmnet",trControl=trainControl(method='repeatedCV',number=5, repeats=10))
  
  vip1<-vip_model(ENET1)
  predxy<-data.frame(x=as.matrix(x), y=y)
  if(!all(is.nan(vip1$Overall))){
    vip1<-  vip1[  vip1$Overall!=0,]
    xnew<-as.matrix(x[,vip1$gene])
    if(ncol(xnew)>10){
      ENET2  <- train(as.matrix(xnew), y=y,
                      method = "glmnet",trControl=trainControl(method='repeatedCV',number=5, repeats=10))
                     # method = "glmnet",tuneGrid = egrid, trControl=trainControl(method='LOOCV'))
    } else {
      ENET2 <- ENET1
    }
    predxy<-data.frame(x=as.matrix(xnew), y=y)
    
    vip2<-vip_model(ENET2)
    vip2<-  vip2[  vip2$Overall>25,] # <------------------------------------------ !!!
    
    bestRes<-ENET2$results[ENET2$results$alpha==ENET2$best$alpha &
                             ENET2$results$lambda==ENET2$best$lambda,]
    
    # final LOO
    predxyLOO<-data.frame()
    for(d in 1:length(y)){
      yloo<-y[-d]
      xloo<-xnew[-d,]
      fit<-glmnet(x=xloo,y=yloo,alpha = bestRes$alpha ,lambda = bestRes$lambda)
      pred<-predict(fit,matrix(xnew[d,], nrow=1))
      predxyLOO<-rbind(predxyLOO,data.frame(d=d,pred=pred,obs=y[d]))
    }
    
    
  } else {
    # in case only intercept estimated, ie no good model can be derived
    bestRes<-data.frame(alpha=NA, lambda=NA, RMSE=NA, Rsquared=NA)
    vip2=NA
  }
  
  return(list(vip=vip2, bestRes=bestRes,ENET2=ENET2,predxy=predxy,predxyLOO=predxyLOO))
}


# Function [wrapper] to do the regression for each of the 4 AEs 
WGCNA.doRegression.allAE<-function(deltacluster=NA, deltagene=NA, reflistClean,  
                                   type,
                                   metric=metric_i){
  
  
  listAE = as.character(unique( reflistClean$ae_descr))
  listScens = expand.grid(listAE,metric)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  out<-apply(listScens,1,function(x,deltacluster,deltagene, reflist){
    ae_i= x[1]
    metric_i = x[2]
    #ae_i="HYPERTROPHY"; metric_i="deltaLFCmean";     #reflist=reflistClean
    if(length(deltagene)>1 ){
      matxy<-mergeObsGeneFC(reflist = reflist, fc_matrix = deltagene, ae=ae_i)
    } 
    if(length(deltacluster)>1 ){
      matxy<-mergeObsClust(reflist = reflist,deltaMetrics = deltacluster,ae=ae_i, metric=metric_i)
    }
    regfit<-fitREG_LOO(x=matxy$x, y=matxy$y)
    names(regfit)<-ae_i
    return(regfit)
  }, reflist=reflistClean,deltacluster=deltacluster,deltagene=deltagene)
  
  stopCluster(cl)
  
  return(list(out=out, scens_ae_metric=listScens))
}



############################################################################
# process Regression Output 


# function to prorcess regression output WGCNA
processClustersRegression<-function(cluster, results.WGCNA){
  
  res<-results.WGCNA[[1]]
  out<-list()
  for(l in 1: length(res)){
    x<-res[[l]][[1]]
    out[[l]]<-list()
    names(out)[l]<-names(res[[l]])[1]
    for(clus in x$gene){
      genes<-as.character(cluster$results$genes[cluster$results$modules == as.character(clus)])
      ngenes<-length(genes)
      nn<-(length(out[[l]])+1)
      out[[l]][[nn]]<-list(genes=genes, ngenes=ngenes)
      names(out[[l]])[nn]<-clus
    }
  }
  return(out)
}

# function to process genes regression output prior to enrichment
cleanGenesOutput<-function(i,d) {
  out<-list(d[[1]][[i]][[1]]$gene)
  names(out)<-names(d[[1]][[i]])[1]
  return(out)
}
# function to process genes regression output prior to enrichment
cleanClustersOutput<-function(i,d) {
  out<-    lapply(d, function(x,i) d[[i]], i=i)
  return(out)
}

#############################################
# enrichment
enrichmentPathwaysSingle<-function(l,genelist,explist,allgenes){
  genelist_i   <- genelist[[l]]
  pathway      <- names(genelist)[l]
  allgeneslist <-length(genelist_i)
  allgene_in_data<-length(explist)
  
  # contingency table
  genes_in_input_list      = length(intersect(explist, genelist_i))
  genes_in_go_term_list    = allgeneslist - genes_in_input_list
  genes_not_in_input_list  = allgene_in_data - genes_in_input_list      
  genes_in_reflist <- allgenes - genes_not_in_input_list -genes_in_go_term_list-genes_in_input_list
  ct   <- matrix(c(genes_in_input_list, genes_in_go_term_list,genes_not_in_input_list ,genes_in_reflist), 
                 byrow = F,ncol=2)
  ft    <- fisher.test(ct,alternative="greater")
  pval  <- ft$p.value
  nlist <- paste(genes_in_input_list,"/",allgeneslist,sep="")
  res   <- data.frame(pathway=pathway, pval=pval, nlist=nlist)
  return(res)
}

# Enrichment for list of different genesets to be enriched (e.g. all significant clusters from one analysis)
enrichmentPathways<-function(genelist, explist, pval=0.01, allgenes=23896){
  allgene_in_data <- length(explist)
  res <-  lapply(1:length(genelist), FUN = enrichmentPathwaysSingle ,
                 genelist=genelist,explist=explist,allgenes=allgenes)
  res<-do.call("rbind", res)
  res<-res[res$pval<pval,]
  res<-res[order(res$pval),]
  return(res)
}

# function to do enrichment given list for different AEs and for each AE 1 or more genesets
doEnrichment<-function(dat.genes, genelist, cluster=T){
  res<-data.frame()
  
  # loop over AE subtypes
  for(l in 1:length(dat.genes)){
    
    ae=names(dat.genes)[[l]]
    
    # loop over clusters (or in case of genes just 1 time)
    for(c in 1:length(dat.genes[[l]])){
      if(cluster) {
        res.enrich<-enrichmentPathways(genelist =genelist, 
                                       explist =  dat.genes[[l]][[c]]$genes) 
        if(nrow(res.enrich)>0){
          res.enrich$ngenes<-unique(dat.genes[[l]][[c]]$ngenes)
          res.enrich$clusname<-unique(names(dat.genes[[l]])[c])
        }
      } else {
        res.enrich<-enrichmentPathways(genelist =genelist, 
                                       explist =  dat.genes[[l]][[c]]) 
        if(nrow(res.enrich)>0){
          res.enrich$ngenes<-NA
          res.enrich$clusname<-NA
        }
      }
      
      if(nrow(res.enrich)>0){
        res.enrich$geneclust=NA
        res.enrich$ae<-ae
        res.enrich$geneclust=c
        res<-rbind(res, res.enrich)
        
      }
      
    }
  }
  return(res)
}

# wrapper function to do the full enrichment
runEnrichment<-function(file.data_processed="rdata/5data_regressionOutput_processed.Rdata",
                        file.reflists="rdata/data_genelists.Rdata",
                        data.genelists=c("genelist.genes.all","genelist.wgcna.all"),
                        ref.genelists=c("genelist.GOBP","genelist.WP","genelist.KEGG","genelist.REACT"),
                        cluster=F ){  # cluster=if genelist true or false
  
  load(file.reflists)
  load(file.data_processed)
  
  # setup list of all combinations of reference lists and gene lists
  enrichlist<-expand.grid(data.genelists=data.genelists, ref.genelists=ref.genelists)
  allres<-data.frame()
  
  for(n in 1:nrow(enrichlist)){
    print(n)
    
    # get pathway list and gene list 
    data.list<-get(as.character(enrichlist$data.genelists[n]))
    reflist<-get(as.character(enrichlist$ref.genelists[n]))
    
    # enrichment; this function loops over multiple clusters and AEs 
    res.enrich<-doEnrichment(dat.genes = data.list, genelist = reflist, cluster=cluster)
    
    # post process
    res.enrich$data.list<-as.character(enrichlist$data.genelists[n])
    res.enrich$gene.list<-as.character(enrichlist$ref.genelists[n])
    res.enrich$data.list<-gsub("genelist.","",res.enrich$data.list)
    res.enrich$method<-sapply(strsplit(res.enrich$data.list,split="\\."),function(x) x[1])
    res.enrich$data<-sapply(strsplit(res.enrich$data.list,split="\\."),function(x) x[2])
    res.enrich$data.list<-NULL
    res.enrich$gene.list<-gsub("genelist.","",res.enrich$gene.list)
    allres<-rbind(allres, res.enrich)
  }
  return(allres)
}



###########################################################
# figure network plots                                    #
###########################################################

# get all clusters for different AEs
getClusterAEs<-function(results.WGCNA){
  res<-lapply(results.WGCNA$out, function(x) {
    ae<-as.character(names(x)[1])
    cluster<-unlist(x[[1]]$gene)
    return(data.frame(ae=ae,cluster=cluster))
  })
  res<-do.call("rbind",res)
  return(res)
}

# function to generate the plots
makeNetworkPlots<-function(cluster, resClusterAE, exclude=NULL,onecluster=NULL){
  #tom<-cluster$TOM
  adj<-cluster$adjacency
  #rownames(tom)<-rownames(adj)
  #colnames(tom)<-colnames(adj)
  
  clusters<-as.character(unique(resClusterAE$cluster))
  count=1
  netplots<-list()
  if(!is.null(onecluster)) {clusters=onecluster}
  for(c in clusters){
    
    whichAE<-paste(resClusterAE$ae[resClusterAE$cluster==c],collapse=", ")
    
    if(!is.null(exclude)){
      aes=resClusterAE$ae[resClusterAE$cluster==c]
      rem<-which(aes==exclude)
      if(length(rem)>0) aes<-aes[-rem]
    }
    if(length(aes)>0){
      
      
      mids<-cluster$results$mID[cluster$results$module==c][1]
      whichMID<-paste("MID",mids)
      title<-paste(whichMID, " (",whichAE,")",sep="")
      
      title<-paste(whichMID, " (",whichAE,")",sep="")
      genes_module<-cluster$results$genes[cluster$results$modules==c]
      #tom1<-tom[genes_module,genes_module]
      #tom1[tom1<median(tom1)]<-0
      
      # ADJ M
      adj1<-adj[genes_module,genes_module]
      adj1[adj1<median(adj1)]<-0
      netplots[[count]]<-ggnet2(adj1, node.size = 6, node.color = "orange",
                                label=TRUE, edge.size = .6, edge.color = "grey")
      netplots[[count]]<-  netplots[[count]]+ggtitle(paste("ADJ",title))
      
      count=count+1
    }
  }
  return(netplots)
}

printNetworkPlots<-function(netplots, type="png", filename="testnetwork"){
  # used to design the plot row x column
  n <- length(netplots)
  nCol <- floor(sqrt(n))
  
  if(type=="png"){png(paste("figures/",filename,".png",sep=""), width=1500, height=1500) }
  if(type=="pdf"){pdf(paste("figures/",filename,".pdf",sep=""), width=18, height=18) }
  do.call("grid.arrange", c(netplots, ncol=nCol))
  dev.off()
}
