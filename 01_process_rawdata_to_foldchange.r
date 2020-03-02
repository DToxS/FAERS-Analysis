rm(list=ls())
##### 
# Fold change data for main drugs set

# Function to convert rawdata files from the raw RNAseq data to be deposited in DTOXS repo to fold change values into a single data.frame# Create list of files to import
files<-dir(path="rawdata/")

# These are drugnames and associated codes used in the RNAseq data
drugNames<-read.csv("reference_files/drugnames_codes.csv", nrow=44, stringsAsFactors = F)
drugNames<-drugNames[grep("nib",drugNames$Drug),]

# remove combos
files<-files[-grep("\\+",files)]

# only keep TKI files
files<-files[grep(paste(drugNames$Code,collapse = "|"),files)]

alldat_fc<-data.frame()
for(f in files){
  print(f)
  file_i <- paste("./rawdata/DEGcor/",f,sep="")
  ff<-read.csv(file_i, sep="\t")
  
  # extract plate and human identifiers
  plate <- substr(unlist(strsplit(f,"Plate."))[[2]],1,1)
  human <- substr(f,7,7)
  code  <- substr(f,35,37)#}
  
  
  # identify drugname from code
  drug  <- drugNames$Drug[drugNames$Code==code]
  
  if(length(drug)>0){
    print(paste(plate,human,code,drug))
    # save to data.frame
    alldat_fc<-rbind(alldat_fc, data.frame(gene   = rownames(ff),
                                           plate  = plate,
                                           human  = human,
                                           drug   = drug,
                                           code   = code, 
                                           logFC  = ff$logFC,
                                           pvalue = ff$PValue,
                                           fdr    = ff$FDR,
                                           logCPM = ff$logCPM))
    
  }
  
}




save(x=alldat_fc, file="rdata/data_rawdata_to_foldchange.Rdata")

#########################################################
# Fold change data for test drugs

rm(list=ls())
 # Function to convert rawdata files from the raw RNAseq data to be deposited in DTOXS repo to fold change values into a single data.frame

 # Create list of files to import
 files<-dir(path="rawdata/Coen-70-Samples-Results/FDR-0.1/")

 # These are drugnames and associated codes used in the RNAseq data
 drugNames<-read.csv("rawdata/drugnames_codes.csv", nrow=44, stringsAsFactors = F)
 drugNames<-drugNames[grep("nib",drugNames$Drug),]



 drugNames<-rbind(drugNames,data.frame(Drug=c("Ibrutinib","Lenvatinib","Nintedanib"),
                                       Code=c("IBR","LEN","NIN"),
                                       Type="TKI"))


 # remove combos
 #files<-files[-grep("\\+",files)]

 # only keep TKI files
 files<-files[grep(paste(drugNames$Code,collapse = "|"),files)]

 alldat_fc2<-data.frame()
 for(f in files){
   print(f)
   file_i <- paste("./rawdata/Coen-70-Samples-Results/FDR-0.1/",f,sep="")
   ff<-read.csv(file_i, sep="\t")

   # extract plate and human identifiers
   plate <- 10
   human <- substr(f,7,7)
   code  <- substr(f,35,37)#}

   
     # identify drugname from code
     drug  <- drugNames$Drug[drugNames$Code==code]

     if(length(drug)>0){
       print(paste(plate,human,code,drug))
       # save to data.frame
       alldat_fc2<-rbind(alldat_fc2, data.frame(gene   = rownames(ff),
                                                plate  = plate,
                                                human  = human,
                                                drug   = drug,
                                                code   = code,
                                                logFC  = ff$logFC,
                                                pvalue = ff$PValue,
                                                fdr    = ff$FDR,
                                                logCPM = ff$logCPM))


   }

 }


 alldat_fc2$plate<-as.factor(alldat_fc2$plate)

 alldat_fc2$drug<-as.character(alldat_fc2$drug)
 alldat_fc2$code<-as.character(alldat_fc2$code)

 alldat_fc2$drug[alldat_fc2$drug=="Ceritinib"]<-"Ceritinib_ext"
 alldat_fc2$code[alldat_fc2$code=="CER"]<-"CER_ext"


 
save(x=alldat_fc2, file="rdata/data_rawdata_to_foldchange_ExtraDrugs.Rdata")
 




