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
```bash
python ./code/parse_gtf.py ./data/'gencode.v12.annotation.gtf' ./output/'gene_annot.parsed.txt'
```
- SNP annotation
First do a quick renaming of the field headers to match the desired column names for downstream processing
```bash
sed -e 's/Chr/chromosome/g' -e 's/Ref_b37/ref_vcf/g' -e 's/Alt/alt_vcf/g' ./data/geuvadis.annot.txt > ./data/snp_annotation.txt
```
>>Split the annotation file into individual chromosomes. You will end up with 22 files for each autosome
```bash
python ./code/split_snp_annot_by_chr.py ./data/snp_annotation.txt ./output/snp_annot
```

- Genotype
```bash
sed 's/Id/varID/g' ./data/geuvadis.snps.txt > ./data/genotype.txt
```
>>Split the genotype into individual chromosomes
```bash
python ./code/split_genotype_by_chr.py ./data/genotype.txt ./output/genotype
```
## Preprocess the expression file
Using your R software load the expression file
```r
# Load packages
library(tidyverse)
library(dplyr)
library(RSQLite)

# Load data
gene_exp = read.table(file = "./data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt", header = TRUE, sep = "\t" )

# Dropping columns we don't need in the gene expression dataframe
gene_exp = gene_exp[,-c(1:4)]

# Name columns (genes)
n = gene_exp$Gene_Symbol
colnames(expr) <- n

# Write the transposed gene expression
write.table(gene_exp_transpose, file = "./output/transformed_expression.txt", sep = "\t",
            row.names = TRUE)
```
### Compute principal components

Make a simple Principal component analysis, and these will be used as covariates to perform multiple linear regressions for each gene in our gene expression matrix and then save the residuals from the regressions as our new expressions for further analysis where needed. This part was formerly [PEER factors](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/), but it has been changed to PCs since it's simpler to interpret, faster to compute and in the end have a good performance, as reported in [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02761-4).

The **PCAForQTL** R library offers two methods to choose the number of PCs to include in the analysis (K): the elbow method and the BE algorithm. According to them "the number of PCs chosen via BE should be considered an upper bound of the reasonable number of PCs to choose in GTEx eQTL data". In this tutorial, the elbow method is used, however, documentation on how to apply the BE algorithm can be found [here](https://github.com/heatherjzhou/PCAForQTL).

```r
library(PCAForQTL)
# Compute principal components
prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE)
PCs<-prcompResult$x # 462 x 462

# Choose K (number of PCs to be used)
## Elbow method
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult) # 29

# Subset PCs, this is the covariates matrix
PCsTop<-PCs[,1:resultRunElbow]

# Save covariates (transposed)
write.table(t(PCsTop), file = "./output/covariates.txt", sep = "\t",
            row.names = TRUE)
```

### Regress out the covariates
You may need the residuals for other analysis e.g estimating h2. We are going to regress out the covariates and save the residuals for any other analysis

```r
# Copy of the expression dataframe to replace values with regressed values
expression = expr

# This loops through all the columns of the transposed gene expression which correspond to each gene, for each gene it runs linear regression on the PCs covariates. Then it sets the residuals to the new expression for that gene.

for (i in 1:length(colnames(expr))) {
  fit = lm(expr[,i] ~ as.matrix(PCsTop))
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
```bash
mkdir -p ./summary ./covariances ./weights
```
- Actual training
Once everything is set you can process a chromosome at a time by executing this code

`Rscript ./code/gtex_tiss_chrom_training.R 1` to train the model for all genes in chromosome 1

To run all the chromosomes you can use a for loop like this
```bash
for i in {1..22}
do
  Rscript ./code/gtex_tiss_chrom_training.R ${i}
done
```
This will take sometime to run and it would be nice to parallelize the step

## Make a database
Make dir for the database
```bash
mkdir -p ./dbs
```
Create database once we have our model summaries we combine them into a single file then create a database
Using R run this chuck of code below
```r
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
```r
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
