#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Compare two gff3 files.
If the gene models of the two gff3 files overlap,
keep the gene model of the first gff3 file.
If they do not overlap, 
keep the gene model of the second gff3 file.
This is beneficial for the screening work of gene family identification.""",
    add_help=False, formatter_class=argparse.RawTextHelpFormatter,
    epilog='Date:2024/09/22 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gff3_file1', metavar='gff3_file1',
                          help='A input file of gff3 format, High precision gff3 file.')
    required.add_argument('gff3_file2', metavar='gff3_file2',
                          help='A input file of gff3 format, Supplementary gff3 file.')
    optional.add_argument('-o', '--output', metavar='str', default=None,  
                          help=f'A output file of gff3 format. defualt=None')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    
    gff3_file1 = args.gff3_file1
    gff3_file2 = args.gff3_file2
    outputfile = args.output
    

def parser_gff3(gff3):
    if gff3.endswith('.gz'):
        f = gzip.open(gff3, 'rt')
    else:
        f = open(gff3, 'r')
        
    genes = {}
    mRNAs = {}
    features = {}
    gene2mRNA = {}

    for l in f:
        if not l.startswith('#') and l.strip() != '':
            l = l.split('\t')
            if l[2] == "gene":
                geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                l[8] = l[8].strip()
                genes.setdefault(geneID, l)

            elif l[2] == "mRNA":
                mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                l[8] = l[8].strip()
                mRNAs.setdefault(mRNAID, l)
                if geneID in gene2mRNA:
                    gene2mRNA[geneID].append(mRNAID)
                else:
                    gene2mRNA.setdefault(geneID, [mRNAID])
            else:
                mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                l[8] = l[8].strip()
                if mRNAID in features:
                    features[mRNAID].append(l)
                else:
                    features.setdefault(mRNAID, [l])
    f.close()
    return genes, mRNAs, features, gene2mRNA

genes1, mRNAs1, features1, gene2mRNA1 = parser_gff3(gff3_file1)
genes2, mRNAs2, features2, gene2mRNA2 = parser_gff3(gff3_file2)

if outputfile == None:
    out = sys.stdout
else:
    out = open(outputfile, 'w')

geneIDs = []
for g1 in genes1:
    for g2 in genes2:
        if genes1[g1][0] == genes2[g2][0] and genes1[g1][6] == genes2[g2][6]:
            if (int(genes1[g1][3]) <= int(genes2[g2][3]) and int(genes1[g1][4]) >= int(genes2[g2][3])) or \
            (int(genes1[g1][3]) >= int(genes2[g2][3]) and int(genes1[g1][3]) <= int(genes2[g2][4])):
                geneIDs.append(g2)
                print('\t'.join(genes1[g1]), file=out)
                for m1 in gene2mRNA1[g1]:
                    print('\t'.join(mRNAs1[m1]), file=out)
                    for f in features1[m1]:
                        print('\t'.join(f), file=out)
                print(file=out)
                
                print('#'+'\t'.join(genes2[g2]), file=out)
                for m2 in gene2mRNA2[g2]:
                    print('#'+'\t'.join(mRNAs2[m2]), file=out)
                    for f in features2[m2]:
                        print('#'+'\t'.join(f), file=out)
                print(file=out)
                
for g2 in genes2:
    if g2 not in geneIDs:
        print('\t'.join(genes2[g2]), file=out)
        for m2 in gene2mRNA2[g2]:
            print('\t'.join(mRNAs2[m2]), file=out)
            for f in features2[m2]:
                print('\t'.join(f), file=out)
        print(file=out)

out.close()
