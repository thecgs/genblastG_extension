#!/usr/bin/env python
# coding: utf-8

import re
import sys
import argparse
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Filter overlap genes.

When two gene models overlap, the model with the highest predicted score is retained.
If the predicted scores of the two gene models are equal, the longest gene model is retained.
""", add_help=False, formatter_class=argparse.RawTextHelpFormatter, 
    epilog='Date:2024/04/23 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('input',  metavar='gff3', help='A input file of gff3 format.')
    optional.add_argument('-o', '--output',  metavar='str', help='A output file of gff3 format.', default=None)
    optional.add_argument('-help', '--help', action='help', help='Show this help message and exit.')
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    gff3_file = args.input
    output = args.output
    
genes = defaultdict(list)
with open(gff3_file) as f:
    for l in f:
        if l.strip() !="" and not l.startswith('#'):
            l = l.strip().split('\t')
            if l[2] == "gene":
                ID = re.search('ID=(.*?);', l[8]).group(1)
                genes[ID] = [l[0], int(l[3]), int(l[4]), l[6], l[5], 0, ID]
            if l[2] == "CDS":
                ID = re.search('Parent=(.*?);', l[8]).group(1).replace('mrna.', '')
                genes[ID][5] += int(l[4]) - int(l[3]) + 1
                
genes = sorted(genes.values(), key=lambda x: (x[0], x[1]))
save_genes = []
for i in range(0, len(genes)):
    if i == 0:
        save_genes.append(genes[i])
    else:
        if (genes[i][0] == save_genes[-1][0]) and (genes[i][3] == save_genes[-1][3]):
            if save_genes[-1][1] <= genes[i][1] <= save_genes[-1][2]:
                if float(genes[i][4]) > float(save_genes[-1][4]):
                    #print(genes[i], save_genes[-1])
                    save_genes[-1] = genes[i]
                elif float(genes[i][4]) == float(save_genes[-1][4]):
                    if float(genes[i][5]) > float(save_genes[-1][5]):
                        #print(genes[i], save_genes[-1])
                        save_genes[-1] = genes[i]
            else:
                save_genes.append(genes[i])
        else:
            save_genes.append(genes[i])

if output == None:
    out = sys.stdout
else:
    out = open(output, 'w')

geneIDs = [gene[6] for gene in save_genes]
with open(gff3_file) as f:
    for l in f:
        if l.strip() !="" and not l.startswith('#'):
            l = l.strip().split('\t')
            if l[2] == "gene":
                ID = re.search('ID=(.*?);', l[8]).group(1)
                if ID in geneIDs:
                    print('\t'.join(l), file=out)
            else :
                ID = re.search('Parent=(.*?);', l[8]).group(1).replace('mrna.', '')
                if ID in geneIDs:
                    print('\t'.join(l), file=out)
out.close()
