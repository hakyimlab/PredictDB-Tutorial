# Turns the parsed genotype annotation file into an RDS file.
# Run from commandline.  First arg is the parsed text file, second is the output file.
argv <- commandArgs(trailingOnly = TRUE)

gene_annot <- read.table(argv[1], stringsAsFactors = FALSE, header = TRUE)
rownames(gene_annot) <- gene_annot$gene_id
saveRDS(gene_annot, argv[2])
