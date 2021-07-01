# PredictDB-Tutorial

## Goals

This tutorial will provide the steps you need to follow to get you started training prediction models and putting them in a format that can be used to run the predixcan suite of software.

## Getting the training dataset
The training dataset can ne downloaded from this [link]().The data set includes the following files you will use for learning;
  1. Normalized gene expression - `GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`
  2. SNP annotation file - `geuvadis.annot.txt`
  3. Gene annotation file - `gencode.v12.annotation.gtf.gz`
  4. Genotype (text format) - `geuvadis.snps.dosage.txt`

Download the files into a directory named `data`.


## Preprocess the genotype
- Parse the gene annotation file
Use the script ` parse_gtf.py` to parse the gene annotation file. This script extracts the important columns about the genes into a file that will be used downstream. The extracted columns include _gene name, gene id, start, end, chr and gene type_.
Run the following command to get the parsed file
```
python ./code/parse_gtf.py ./data/'gencode.v12.annotation.gtf' ./output/'gene_annot.parsed.txt'
```
- SNP annotation
First do a quick renaming of the field headers to match the desired column names for downstream processing
```
sed -e 's/Chr/chromosome/g' -e 's/Ref_b37/ref_vcf/g' -e 's/Alt/alt_vcf/g' ./data/geuvadis.annot.txt > ./data/snp_annotation.txt
```
Split the annotation file into individual chromosomes. You will end up with 22 files for each autosome
```
python ./code/split_snp_annot_by_chr.py ./data/geuvadis.annot.txt ./output/snp_annot
```

## Preprocess the expression file
Using your R software load the expression file
```
# Load packages
library(tidyverse)
library(dplyr)
library(RSQLite)

# Load data
gene_exp = read.table(file = "./data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt", header = TRUE, sep = "\t" )

##Dropping columns we don't need in the gene expression dataframe
gene_exp = gene_exp[-c(1, 3, 4)]

# rename a column
gene_exp = rename(gene_exp, 'Gene_Name' = Gene_Symbol)

# Transpose the data to enable you rub peer factors
n = gene_exp$Gene_Name
gene_exp_transpose <- as.data.frame(t(gene_exp[,-1]))
colnames(gene_exp_transpose) <- n

# write the transposed gene expression. The .csv is for the peer tool
write.table(gene_exp_transpose, file = './output/gene_exp.csv', sep = ",", col.names = TRUE, row.names = FALSE)
write.table(gene_exp_transpose, file = "./output/transformed_expression.txt", sep = "\t",
            row.names = TRUE)
```
- Calculate peer covariates

Generate [PEER factors](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/), these peer factors will be used as covariates to perform multiple 
linear regression for each gene in our gene expression matrix and then save the residuals from the regression as our new expressions for further analysis where needed.

Ensure you have `peer` tool installed on your machine. There is a description of how to download the PEER tool [here](https://github.com/hakyimlab/peer).

Use the PEER tool to generate PEER factors from our transformed gene expression file (`./output/gene_exp.csv`). According to GTEx protocol, If the number of samples is greater than or equal to 350, we use 60 PEER factors. If the number of samples is between 250 and 350, we use 45. Between 150 and 250, we use 30, and less than 150 we use 15. For this study, the number of samples is 463 so we will use 60 PEER factors.

```
peertool -f './output/gene_exp.csv' -n 60 --has_header -o ./output/peer_out
```
__Note:__ this takes a long time to run the peers

Once completed read the output of peertool into r to add column names. This will form our covariates matrix
```
peer_factors = read.csv(file = "./output/peer_out/X.csv", header = FALSE)

#Set the column names for the PEER factors (covariates) as the subject IDs
colnames(peer_factors) = rownames(gene_exp_transpose)

# write out a covariates matrix
write.table(peer_factors, file = "./output/covariates.txt", sep = "\t",
            row.names = TRUE)
```
- Regress out the covariates
You may need the residuals for other analysis e.g estimating h2. We are going to regress out the covariates and save the residuals for any other analysis

```
## Make a copy of the transposed gene expression dataframe so that we can replace the values with the residuals of the multiple linear regressions.
expression = gene_exp_transpose

# This loops through all the columns of the transposed gene expression which correspond to each gene,
for each gene it runs linear regression on the PEER factor covariates. Then it sets the residuals to the new expression for that gene.

for (i in 1:length(colnames(gene_exp_transpose))) {
    fit = lm(gene_exp_transpose[,i] ~ t(as.matrix(peer_factors)))
    expression[,i] <- fit$residuals
  }

# Write out the residual expression file
write.table(expression, file = "./output/residuals_expression.txt", sep = "\t",
            row.names = TRUE)
```
## Train the models


or follow this nextflow pipeline



## Data
Data for this tutorial is from [here](https://uchicago.box.com/s/ewnrqs31ivobz2sn6462cq2eb423dvpr)

- [ ] TODO: write a post in lab-notes about how to install nextflow and a brief description on how to use it. More details specific to cri if needed.




## This is a tutorial for using GEUVADIS data to run the PredictDB V7 Pipeline
  A description of the pipeline can be found here https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7
