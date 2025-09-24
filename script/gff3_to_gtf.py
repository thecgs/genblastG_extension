#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse

def gff2gtf(gff3_file, gtf_file):
    
    if gff3_file.endswith('.gz'):
        gff3 = gzip.open(gff3_file, 'rt')
    elif gff3_file == '-':
        gff3 = sys.stdin
    else:
        gff3 = open(gff3_file, 'r')
    
    if gtf_file.endswith('.gz'):
        gtf = gzip.open(gtf_file, 'wt')
    elif gtf_file == '-':
        gtf = sys.stdout
    else:
        gtf = open(gtf_file, 'w')
    
    gene = None
    count = 0
    for l in gff3:
        if (not l.startswith('#')) and (l.strip() != ''):
            l = l.strip()
            if l[-1] != ';':
                l += ";"
            l = l.split('\t')
            
            if l[2] == 'gene':
                gene = "\t".join(l[:8])
                
            elif l[2] == 'mRNA':
                if count !=0:
                    print(file=gtf)
                count += 1
                l[2] = 'transcript'
                gene_id = re.search(f"Parent=(.*?);", l[8]).group(1)
                transcript_id = re.search(f"ID=(.*?);", l[8]).group(1)
                l[8] = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
                
                if gene !=None:
                    print(gene + '\t' + f'transcript_id "{transcript_id}"; gene_id "{gene_id}";', file=gtf)
                print('\t'.join(l), file=gtf)
                
            elif l[2] == 'three_prime_UTR':
                l[2] = '3UTR'
                l[8] = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
                print('\t'.join(l), file=gtf)
                
            elif l[2] == 'five_prime_UTR':
                l[2] = '5UTR'
                l[8] = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
                print('\t'.join(l), file=gtf)
            else:
                l[8] = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
                print('\t'.join(l), file=gtf)
                
    gff3.close()
    gtf.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gff3 format convert to gtf format.', add_help=False, 
                                     epilog='Date:2025/09/24 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gff3', metavar='gff3', help='A file of gff3 format.')
    required.add_argument('gtf', metavar='gtf', help='A file of gtf format')
    optional.add_argument('-h','--help',action='help', help='Show this help message and exit')
    optional.add_argument('-v','--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    gff2gtf(gff3_file=args.gff3, gtf_file=args.gtf)
