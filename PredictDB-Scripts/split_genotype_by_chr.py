#! /usr/bin/env python

import os
import sys

'''
Script to split a GTEx genotype file into multiple files by chromosome.

From commandline, first argument is the genotype file, second is the
prefix for the output files.  The suffix 'chrN.txt' will be added to the
prefix provided, where N is the chromosome number.

In splitting, script will only keep unambiguously stranded SNPs. I.e.,
no INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G
and vice-versa.

The input file is expected to be a tab-delimited text file including a
header row, where the header field is ID (for snp varID) and then a
variable number of fields with the sample ID numbers.  The first column
has the snp ID in the format of (chr_pos_refAll_effAll_build) and the
dosages are encoded on a 0-2 scale representing the number or imputed
number of the effect alleles the sample posseses.
'''

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def split_genotype(geno_file, out_prefix):
    # Make output file names from prefix.
    geno_by_chr_fns = [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]
    # Open connection to each output file.
    geno_by_chr = [open(f, 'w') for f in geno_by_chr_fns]

    with open(geno_file, 'r') as geno:
        # Write header in each file
        header = geno.readline()
        snps = set()
        for f in geno_by_chr:
            f.write(header)

        for line in geno:
            # First attribute of line is is chr_pos_refAllele_effAllele_build
            # Extract this attribute and parse into list
            varID_list = (line.split()[0].split('_'))
            chr = varID_list[0]
            refAllele = varID_list[2]
            effectAllele = varID_list[3]
            # Skip non_single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                continue
            varID = '_'.join(varID_list)
            # Some snps have 2 rows for some reason. Attributes are nearly
            # identical. Only keep the first one found.
            if varID in snps:
                continue
            snps.add(varID)
            # Write line to appropriate file
            index = int(chr) - 1
            geno_by_chr[index].write(line)

    for f in geno_by_chr:
        f.close()

if __name__ == '__main__':
    genotype_file = sys.argv[1]
    out_prefix = sys.argv[2]
    split_genotype(genotype_file, out_prefix)
