# FAERS-Analysis
Source codes for "Transcriptomic profiling of human cardiac cells predicts protein kinase inhibitor-associated cardiotoxicity" 
 
## Script for generation of FAERS risk score signature

### All software dependencies and operating systems (including version numbers)
platform       x86_64-pc-linux-gnu         
arch           x86_64                      
os             linux-gnu                   
system         x86_64, linux-gnu           
status                                     
major          3                           
minor          6.1                         
year           2019                        
month          07                          
day            05                          
svn rev        76782                       
language       R                           
version.string R version 3.6.1 (2019-07-05)
nickname       Action of the Toes   

R package versions used:
caret      6.0-84
glmnet      2.0-18

###Installatin guide

1. Download all raw data from GEO repository
2. Run R and R packages according to listed versions
3. Run R code in sequence the files have been numbered.



### Overview of R scripts:
00_functions.R: contains all functions used
01_process_rawdata_to_foldchange.R: processing of raw data to fold change dataset used for elastic net analysis
02_analysis_riskscores.R: script to compute FAERS risk scores
03_analysis_elastic_net.R: script to perform elastic net regression analysis
04_analysis_test_drugs.R: script to generate predictions for test drugs

### Overview of reference files:
aersMineExploreDataSet_1957.TSV: Raw FAERS data download from AersMine
aersMineExploreDataSet_1957.aers: Description/time stamp of FAERS data download from AersMine
drugnames_codes.xlsx: Listing of drug codes and drug names
meddra_annotate.xlsx: MEDDRA adverse event ontology used for FAERS risk score calculation
vigiAccessTKI.xlsx: ViciAccess Database download for risk score analysis




