#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import gzip
import shutil
import subprocess
import argparse
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""This script is mainly the wrapper of genblastg software.
It can run genblastg in multiple threads, reconstruct the results, and output the standard gff3 file. 
Secondly, it can filter the redundant gene model according to the prediction score and the length of 
the predicted gene to generate the best non redundant gene model.""", 
    add_help=False, formatter_class=argparse.RawTextHelpFormatter,
    epilog='Date:2025/09/24 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-g','--genome', metavar='str',
                          help='A file of genome fasta format.', required=True)
    required.add_argument('-q','--query', metavar='str',
                          help='A file of query protein fasta format.',required=True)
    optional.add_argument('-p','--prefix', metavar='str',
                          help='A prefix of output. default=genblastG_extension', default='genblastG_extension')
    optional.add_argument('-c', '--query_cover', metavar='float', default=0.8, type=float,
                          help='minimum query cover (0-1) to report an alignment. defualt=0.8')
    optional.add_argument('-gap', '--gap', action='store_true', 
                          help="parameter for blast: Perform gapped alignment. default=False")
    optional.add_argument('-e', '--evalue', metavar='str', 
                          help='A maximum evalue of report alignments. defualt=1e-5', default='1e-5')
    optional.add_argument('-G', '--genetic_code', metavar='int', type=int, default=1, 
                          help="Genetic code. default=1")
    optional.add_argument('-t', '--thread', metavar='int', default=16, type=int, 
                          help='Thread number of single sortware. defualt=16')
    optional.add_argument('-h','--help', action='help',
                          help='Show this help message and exit.')
    optional.add_argument('-v','--version', action='version', version='v1.00', 
                          help="Show program's version number and exit.")
    args = parser.parse_args()

def rm(file):
    if os.path.isdir(file):
        shutil.rmtree(file, ignore_errors=True)
    else:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass
    return None

def mkdir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    return None

def split_seq(inputfile, outdir="00.split_seqs"):
    if inputfile.endswith('.gz'):
        handle = gzip.open(inputfile, 'rt')
    elif inputfile == '-':
        handle = sys.stdin
    else:
        handle = open(inputfile, 'r')
    
    mkdir(outdir)
    os.chdir(outdir)
    
    inputfiles = []
    for n, record in  enumerate(SeqIO.parse(handle, 'fasta')):
        with open(f'{record.id}.{str(n+1)}.fasta', 'w') as f:
            seqence = ''
            for i in record.seq:
                if i.islower():
                    seqence += 'N'
                else:
                    seqence += i
                    
            print(f'>{record.id}.{str(n+1)}\n{seqence}', file=f)
            inputfiles.append(os.path.realpath(f'{record.id}.{str(n+1)}.fasta'))
    handle.close()
    os.chdir('../')
    return inputfiles

def run_gblast(inputfile, genome_file, evalue='1e-10', query_cover=0.8, gap=False, output_file='gene_gff', thread=10):
    if gap == False:
        gap = "F"
    else:
        gap = "T"

    alignscore = open('alignscore.txt','a')
    score = open(f"{sys.path[0]}/bin/alignscore.txt",'r')
    alignscore.writelines(score)
    alignscore.close()
    score.close()
    mkdir("01.gblastG_result")
    
    def run(args):
        if os.path.exists('01.gblastG_result/' + args[1] + '_1.1c_2.3_s1_0_16_1.gff'):
            pass
        else:
            cline = subprocess.Popen(args[0], env={'GBLAST_PATH':f"{sys.path[0]}/bin", 'PATH':f"{sys.path[0]}/bin"}, shell=True)
            try:
                cline.communicate(timeout=300)
            except:
                print(f'error:{args[1]}')
            ret_code = cline.poll()
            if ret_code !=0:
                cline.terminate()
        return None
    
    files = split_seq(inputfile, outdir="00.split_seqs")
    args = []
    for file in files:
        cmd = f'genblastG -P blast -q {file} -t {genome_file} -e {evalue} -c {query_cover} -g {gap} -r 5 -gff -o ./01.gblastG_result/{os.path.splitext(os.path.basename(file))[0]}'
        arg = (cmd, os.path.splitext(os.path.basename(file))[0], os.path.basename(genome_file))
        args.append(arg)
        
    with ThreadPoolExecutor(max_workers=thread) as pool:
        pool.map(run, args)

    out = open(output_file, 'w')
    for file in files:
        with open("01.gblastG_result/" + os.path.splitext(os.path.basename(file))[0] + '_1.1c_2.3_s1_0_16_1.gff', 'r') as f:
            out.write(f.read())
    out.close()
    return None

def get_phases(pos):
    phases = [0]
    for i, p in enumerate(pos[1:]):
        phase=(3 - ((pos[i][1] - pos[i][0] + 1) - phases[i]) % 3) % 3
        phases.append(phase)
    return phases

def reconsitution_gff(inputfile, outputfile):
    genes = defaultdict(str)
    mRNAs = defaultdict(str)
    CDSs = defaultdict(list)
    count = 1
    with open(inputfile, 'r') as f:
        for l in f:
            if (not l.startswith('#')) and (l.strip() != ''):
                l = l.strip().split('\t')
                if l[2] == 'transcript':
                    GeneID = "GENE" + str(count) + "_" + re.search("ID=(.*?);", l[8]).group(1)
                    count += 1
                    Target = re.search("Name=(.*)", l[8]).group(1)
                    gene = (l[0], l[1], 'gene', l[3], l[4], l[5], l[6], l[7], f"ID={GeneID};Target={Target};")
                    mRNA = (l[0], l[1], 'mRNA', l[3], l[4], '.', l[6], l[7], f"ID=mrna.{GeneID};Parent={GeneID};")
                    genes[GeneID] = '\t'.join(gene)
                    mRNAs[GeneID] = '\t'.join(mRNA)
                elif l[2] == 'coding_exon':
                    CDS = [l[0], l[1], 'CDS', l[3], l[4], '.', l[6], l[7], f"ID=cds.{GeneID};Parent=mrna.{GeneID};"]
                    CDSs[GeneID].append(CDS)
    out = open(outputfile, 'w')
    for ID in genes:
        print(genes[ID], file=out)
        print(mRNAs[ID], file=out)
        
        if mRNAs[ID].split('\t')[6] == '+':        
            CDSs[ID] = sorted(CDSs[ID], key=lambda x: int(x[3]), reverse=False)
            pos = [(int(e[3]), int(e[4])) for e in CDSs[ID]]
            phases = get_phases(pos)
        else:
            CDSs[ID] = sorted(CDSs[ID], key=lambda x: int(x[3]), reverse=True)
            pos = [(int(e[3]), int(e[4])) for e in CDSs[ID]]
            phases = get_phases(pos)
            
        for i, e in enumerate(CDSs[ID]):
            print(e[0], e[1], 'exon', e[3], e[4], e[5], e[6], e[7], f"ID=exon.{ID};Parent=mrna.{ID};",sep='\t', file=out)
            e[7] = str(phases[i])
            print('\t'.join(e), file=out)
        print(file=out)
    out.close()
    return None
run_gblast(inputfile=args.query, genome_file=args.genome,
           evalue=args.evalue, query_cover=args.query_cover, gap=args.gap,
           output_file="genblastG.gff3", thread=args.thread)
reconsitution_gff(inputfile="genblastG.gff3", outputfile=args.prefix+'.raw.gff3')

rm('00.split_seqs')
rm('01.gblastG_result')
rm('genblastG.gff3')
rm('alignscore.txt')
rm('formatdb.log')
rm('perform.txt')

subprocess.run(f"{os.path.join(sys.path[0], 'script/filter_overlap_gene.py')} {args.prefix}.raw.gff3 | {os.path.join(sys.path[0], 'script/sort_gff3.py')} - -o {args.prefix}.filtered.gff3", shell=True)
subprocess.run(f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {args.prefix}.filtered.gff3 {args.genome} -t CDS -o {args.prefix}.cds.fasta -G {args.genetic_code}", shell=True)
subprocess.run(f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {args.prefix}.filtered.gff3 {args.genome} -t prot -o {args.prefix}.pep.fasta -G {args.genetic_code}", shell=True)
