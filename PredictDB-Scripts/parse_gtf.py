#! /usr/bin/env python

import sys

'''
Parses a gtf file to pull chromosome number, gene_id, gene_name, gene start
position, and gene end position, and put in a tab delimited file.

Runs from command line.  First argument is the gtf file, second is the
output file.
'''

HEADER_FIELDS = ['chr', 'gene_id', 'gene_name', 'start', 'end', 'gene_type']

def parse_gtf(gtf_path, out_path):
    with open(gtf_path, 'r') as gtf, open(out_path, 'w') as out:
        out_header = '\t'.join(HEADER_FIELDS) + '\n'
        out.write(out_header)
        for line in gtf:
            # Skip comments.
            if line[0] == '#':
                continue
            gene_fields = line.split('\t')
            # Exclude any rows not describing genes.
            if gene_fields[2] != 'gene':
                continue
            gene_attributes = gene_fields[-1].split('; ')
            attr_dict = dict(attribute.split() for attribute in gene_attributes if attribute)
            # Some gtf files may have chromosome number with 'chr' prefix.
            # We just need the number.
            chr = gene_fields[0][3:] if "chr" in gene_fields[0] else gene_fields[0]
            # Extract rest of desired fields.
            start = gene_fields[3]
            end = gene_fields[4]
            id = attr_dict['gene_id'].strip('"')
            name = attr_dict['gene_name'].strip('"')
            type = attr_dict['gene_type'].strip('"')
            # Concatenate together and write out to file.
            out_line = '\t'.join([chr, id, name, start, end, type]) + '\n'
            out.write(out_line)

if __name__ == '__main__':
    gtf_path = sys.argv[1]
    out_path = sys.argv[2]
    parse_gtf(gtf_path, out_path)
