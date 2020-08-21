source("../code/gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]

#tiss <- argv[1]
#chrom <- argv[2]

snp_annot_file <- "../output/snp_annot.chr" %&% chrom %&% ".txt"
gene_annot_file <- "../output/gene_annot.parsed.txt"
genotype_file <- "../output/genotype.chr" %&% chrom %&% ".txt"
expression_file <- "../output/transformed_expression.txt"
covariates_file <- "../output/covariates.txt"
prefix <- "Model_training"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


