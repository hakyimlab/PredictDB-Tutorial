# Turns a set of 22 snp annotation files (organized by chromosome number)
# into RDS files for fast reading into R.  Runs from the commandline, and
# the first argument is the prefix for the input files.  The prefix provided
# should include everything up to the chromosome number. The output file name
# will be the prefix with the chromosome number and '.RDS' concatenated at
# the end.

"%&%" <- function(a,b) paste(a, b, sep = "")
argv <- commandArgs(trailingOnly = TRUE)

file_head <- argv[1]

for (i in 1:22) {
  infile <- file_head %&% i %&% ".txt"
  outfile <- file_head %&% i %&% ".RDS"
  snp_annot <- read.table(infile, stringsAsFactors = FALSE, header = TRUE)
  snp_annot <- unique(snp_annot)
  rownames(snp_annot) <- snp_annot$varID
  saveRDS(snp_annot, outfile)
  print(i %&% "/22")
}
