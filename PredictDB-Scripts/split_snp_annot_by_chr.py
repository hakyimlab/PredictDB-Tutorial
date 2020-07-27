#! /usr/bin/env python

import sys

'''
Script to split a SNP annotation file into multiple files by chromosome.

From commandline, first argument is the snp annotation file, second is
the prefix for the output files.  The suffix 'chrN.txt' will be added
to the prefix provided, where N is the chromosome number.

In splitting, script will only keep unambiguously stranded SNPs. I.e.,
no INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G
and vice-versa.

Input snp annotation file is expected to be a tab-delimited text file,
with a header row, with fields for chromosome, position, variantID,
reference allele, alternative allele, rsid_label1, rsid_label2, and 
number of alternative alleles per site. See file
GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz
from gtexportal.org for an example of such a file.

The output files will be tab-delimited text files with chromosome,
position, variantID, reference allele, effect allele, and rsid.
NOTE: the rsid number chosen is from rsidlabel2.
'''

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
HEADER_FIELDS = ['chr','pos','varID','refAllele','effectAllele','rsid']

def split_snp_annot(annot_file, out_prefix):
    # Make output file names from prefix.
    snps_by_chr_files= [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]
    # Open connection to each output file
    snp_by_chr = [open(f, 'w') for f in snps_by_chr_files]
    # Write header in each file.
    header = '\t'.join(HEADER_FIELDS)+'\n'
    for f in snp_by_chr:
        f.write(header)
    with open(annot_file, 'r') as ann:
        # Skip header from input file
        ann.readline()
        # Extract rows from input and write to body in appropriate output.
        for line in ann:
            attrs = line.split()
            chr = attrs[0]
            pos = attrs[1]
            varID = attrs[2]
            refAllele = attrs[3]
            effectAllele = attrs[4]
            rsid = attrs[6]
            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                continue
            if rsid == '.':
                continue
            index = int(chr) - 1
            row = '\t'.join([chr,pos,varID,refAllele,effectAllele,rsid])+'\n'
            snp_by_chr[index].write(row)
    # Close connection to each output file.
    for f in snp_by_chr:
        f.close()

if __name__ == '__main__':
    annot_file = sys.argv[1]
    out_prefix = sys.argv[2]
    split_snp_annot(annot_file, out_prefix)
