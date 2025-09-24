#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse
from collections import defaultdict

def gtf2gff(gtf_file, gff3_file):
    
    if gtf_file.endswith('.gz'):
        gtf = gzip.open(gtf_file, 'rt')
    elif gtf_file == '-':
        gtf = sys.stdin
    else:
        gtf = open(gtf_file, 'r')
        
    if gff3_file.endswith('.gz'):
        gff3 = gzip.open(gff3_file, 'wt')
    elif gff3_file == '-':
        gff3 = sys.stdout
    else:
        gff3 = open(gff3_file, 'w')
    
    genes = defaultdict(str)
    mRNAs = defaultdict(str)
    genes2mRNAs = defaultdict(list)
    features = defaultdict(list)
    
    for l in gtf:
        if (not l.startswith('#')) and (l.strip() != ''):
            l = l.strip()
            if l[-1] != ';':
                l += ";"
            l = l.split('\t')
            
            gene_id = re.search('gene_id "(.*?)"', l[8]).group(1)
            transcript_id = re.search(f'transcript_id "(.*?)"', l[8]).group(1)
            
            if l[2] == 'gene':
                genes[gene_id] =  "\t".join(l[:8]) + f"\tID={gene_id};"
                
            elif l[2] == 'transcript':
                l[2] = 'mRNA'
                l[8] = f'ID={transcript_id};Parent={gene_id};'
                mRNAs[transcript_id] = "\t".join(l)
                genes2mRNAs[gene_id].append(transcript_id)
                
            elif l[2] == '3UTR':
                l[2] = 'three_prime_UTR'
                l[8] = f'ID={transcript_id}.utr3;Parent={transcript_id};'
                features[transcript_id].append("\t".join(l))

            elif l[2] == '5UTR':
                l[2] = 'five_prime_UTR'
                l[8] = f'ID={transcript_id}.utr5;Parent={transcript_id};'
                features[transcript_id].append("\t".join(l))
                
            elif l[2] == 'exon':
                l[8] = f'ID={transcript_id}.exon;Parent={transcript_id};'
                features[transcript_id].append("\t".join(l))
                
            elif l[2] == 'CDS':
                l[8] = f'ID={transcript_id}.cds;Parent={transcript_id};'
                features[transcript_id].append("\t".join(l))
                
            else:
                l[8] = f'ID={transcript_id}.{l[2]};Parent={transcript_id};'
                features[transcript_id].append("\t".join(l))
    gtf.close()
    
    if len(genes) != 0:
        for g in genes:
            print(genes[g], file=gff3)
            for m in genes2mRNAs[g]:
                print(mRNAs[m], file=gff3)
                for f in features[m]:
                    print(f, file=gff3)
            print(file=gff3)
    else:
        for g in genes2mRNAs:
            start = min([int(mRNAs[k].split('\t')[3]) for k in genes2mRNAs[g]])
            end = max([int(mRNAs[k].split('\t')[4]) for k in genes2mRNAs[g]])
            gene = mRNAs[genes2mRNAs[g][0]].split('\t')
            gene[3] = str(start)
            gene[4] = str(end)
            gene[2] = 'gene'
            gene[8] = f"ID={g};"
            print('\t'.join(gene), file=gff3)
            for m in genes2mRNAs[g]:
                print(mRNAs[m], file=gff3)
                for f in features[m]:
                    print(f, file=gff3)
            print(file=gff3)
    gff3.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gtf format convert to gff3 format.', add_help=False, 
                                     epilog='Date:2025/09/24 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gtf', metavar='gtf', help='A file of gtf format')
    required.add_argument('gff3', metavar='gff3', help='A file of gff3 format.')
    optional.add_argument('-h','--help',action='help', help='Show this help message and exit')
    optional.add_argument('-v','--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    gtf2gff(gtf_file=args.gtf, gff3_file=args.gff3)
