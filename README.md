# PredictDB-Tutorial

## Goals

This tutorial will provide the steps you need to follow to get you started training prediction models and putting them in a format that can be used to run the predixcan suite of software.

## Getting the training dataset
The training dataset can ne downloaded from this [link](https://uchicago.box.com/s/ewnrqs31ivobz2sn6462cq2eb423dvpr).The data set includes the following files you will use for learning;
  1. Normalized gene expression - `GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`
  2. SNP annotation file - `geuvadis.annot.txt`
  3. Gene annotation file - `gencode.v12.annotation.gtf.gz`
  4. Genotype (text format) - `geuvadis.snps.dosage.txt`

Download the files into a directory named `data` and decompress them.


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
>>Split the annotation file into individual chromosomes. You will end up with 22 files for each autosome
```
python ./code/split_snp_annot_by_chr.py ./data/snp_annotation.txt ./output/snp_annot
```

- Genotype
```
sed 's/Id/varID/g' ./data/geuvadis.snps.txt > ./data/genotype.txt
```
>>Split the genitype into individual chromosomes
```
python ./code/split_genotype_by_chr.py ./data/genotype.txt ./output/genotype
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
- Check all processed data is available  in the output dir
To train the model you need the files you generated from the steps above;
  - gene expression file: `transformed_expression.txt`
  - covariates file: `covariates.txt`
  - residual file after regressing out covariates: `residuals_expression.txt`
  - gene annotation file: `gene_annot.parsed.txt`
  - snp_annotation file: `snp_annot.chr*.txt` for chr 1-22
  - genotype file: `genotype.chr*.txt` for chr 1-22

Confirm you have all these files in your output directory.

- Set up the execution script
We will process each chromosmome individually using this rscript `.code/gtex_tiss_chrom_training.R`.
Edit the file to fit the paths to the constant files and the ones which will be changed for each chromosome i.e snp_annotation and genotype files.

- Create the required directories
```
mkdir -p ./summary ./covariances ./weights
```
- Actual training
Once everything is set you can process a chromosome at a time by executing this code

`Rscript ./code/gtex_tiss_chrom_training.R 1` to train the model for all genes in chromosome 1

To run all the chromosomes you can use a for loop like this
```
for i in {1..22}
do
  Rscript ./code/gtex_tiss_chrom_training.R ${i}
done
```
This will take sometime to run and it would be nice to parallelize the step

## Make a database
Make dir for the database
```{bash}
mkdir -p ./dbs
```
Create database once we have our model summaries we combine them into a single file then create a database
Using R run this chuck of code below
```{r}
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
model_summaries <- read.table('./summary/Model_training_chr1_model_summaries.txt',header = T, stringsAsFactors = F)
tiss_summary <- read.table('./summary/Model_training_chr1_summary.txt', header = T, stringsAsFactors = F)
  
n_samples <- tiss_summary$n_samples
  
for (i in 2:22) {
  model_summaries <- rbind(model_summaries,
                            read.table('./summary/Model_training_chr'%&%as.character(i) %&% '_model_summaries.txt', header = T, stringsAsFactors = F))
  tiss_summary <- rbind(tiss_summary,
                             read.table('./summary/Model_training_chr' %&% as.character(i) %&% '_summary.txt', header = T, stringsAsFactors = F))
  
}
  
model_summaries <- rename(model_summaries, gene = gene_id)

# Create a database connection
conn <- dbConnect(drv = driver, './dbs/gtex_v7_models.db')
dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

# Weights Table -----
weights <- read.table('./weights/Model_training_chr1_weights.txt', header = T,stringsAsFactors = F)
for (i in 2:22) {
  weights <- rbind(weights,
              read.table('./weights/Model_training_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
  
}
  
weights <- rename(weights, gene = gene_id)
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = 'EUR') # Provide the population info
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  
# Construction Table ----
construction <- tiss_summary %>%
                    select(chrom, cv_seed) %>%
                    rename(chromosome = chrom)

dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
dbDisconnect(conn)
```

## Filter the database
Filter the database to select significant models using R
```{r}
unfiltered_db <- './dbs/gtex_v7_models.db'
filtered_db <- './dbs/gtex_v7_models_filtered_signif.db'
driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, unfiltered_db)
out_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.1')
model_summaries <- model_summaries %>% 
                    rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
model_summaries$pred.perf.qval <- NA
dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
construction <- dbGetQuery(in_conn, 'select * from construction')
dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)
weights <- dbGetQuery(in_conn, 'select * from weights')
weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
dbDisconnect(in_conn)
dbDisconnect(out_conn)
```
Your database is ready to use.


**_This pipeline has also been implemented in nextflow pipeline [here](https://github.com/hakyimlab/PredictDb-nextflow)_**


- [ ] TODO: write a post in lab-notes about how to install nextflow and a brief description on how to use it. More details specific to cri if needed.




## This is a tutorial for using GEUVADIS data to run the PredictDB V7 Pipeline
  A description of the pipeline can be found here https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7
