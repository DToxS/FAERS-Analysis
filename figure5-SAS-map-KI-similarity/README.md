# Code for Figure 5 

## Kinase Selectivity figures 

`kinase-selectivity.R` contains the code used to generate subfigures seen in Figure 5. 

To run, first install dependancies: 

`source("https://bioconductor.org/biocLite.R")`

`biocLite("BiocInstaller")`

`install.packages(c("data.table", "devtools", "magrittr" , "ggplot2", "scales" , "tidyr" , "ggrepel", "caret", "glmnet" ))`

`devtools::install_github("ropensci/iheatmapr")`

Then run the R script interactively in R studio/VSCode 

## Generation of Tanimoto coefficients 

`pairwise_ecfp4_distances.py` contains code used to generate pairwise tanimoto coefficients (Tc) for the kinase inhibitors tested 

To run, first install the RDKit dependancy: 

Install miniconda: https://docs.conda.io/en/latest/miniconda.html 

Follow installation procedure to install RDkit: https://www.rdkit.org/docs/Install.html

To run: 

`python3 pairwise_ecfp4_distances.py [SMILES FILE] > pairwise-distances.csv`

Example run: 

`python3 pairwise_ecfp4_distances.py data/KI.smi > pairwise-distances.csv` 

The script will generate a table in comma separated format called `pairwise-distances.csv` that will have the following: 

Column1 D1 (Drug 1)

Column2 D2 (Drug 2)

Column3 ECFP4 (Drug1/2 ECFP radius 2 fingerprint Tc)

Column4 ECFP2 (Drug1/2 ECFP radius 1 fingerprint Tc)

Column5 MACCS (Drug1/2 MACCs fingerprint Tc)

Column6 DL (Daylight fingerprint Tc)

Column7 AVG (Average fingerprint Tc)

Column8 Weighted (Weighted average fingerprint Tc )

*Weights are as follows: 30% ECFP4 + 30%ECFP2 + 30%DL + 10%MACCS* 







