#!/usr/bin/env python
# coding: utf-8


import os
import sys
import csv
import subprocess
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='By genblastG software, homologous gene of target species was predicted based on homologous protein sequences, and GFF files was generated.',usage='python run_genblastG.py -i seq.fa -g genome.fa -o out.gff3',add_help=False,epilog='date:2020/08/26 author:guisen chen email:thecgs001@foxmail.com')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-g','--genome',metavar='[genome.fa]',help='genome.fa',required=True)
required.add_argument('-i','--input',metavar='[seq.fa]',help='seq.fa',required=True)
required.add_argument('-o','--output',metavar='[out.gff3]',help='out.gff3',required=True)
optional.add_argument('-h','--help',action='help',help='show this help message and exit')
optional.add_argument('-v','--version',action='version',version='v1.00')
args = parser.parse_args()

def get_single_seq(seq_file):
    seq_records = open(seq_file,'r')
    id_list = []
    for seq_record in SeqIO.parse(seq_records,'fasta'):
        id_list.append(seq_record.id)
        single_seq_file = open(f'{seq_record.id}.fa','w')
        single_seq_file.write(f'>{seq_record.id}\n{seq_record.seq}')
        single_seq_file.close()
    seq_records.close()
    return id_list
    
def run_gblast(single_seq,genome_file):
    cline = subprocess.Popen(f'genblastG -P blast -q {single_seq} -t {genome_file} -e 1e-10 -r 5 -c 0.8 -gff -o gblast.out',env={'GBLAST_PATH':f"{sys.path[0]}",'PATH':f"{sys.path[0]}"},shell=True)
    try:
        cline.communicate(timeout=300)
    except:
        print(f'error:{single_seq}')
    ret_code = cline.poll()
    if ret_code !=0:
        cline.terminate()
    gene_gff = open('gene_gff','a')
    headle_gff = open('gblast.out_1.1c_2.3_s1_0_16_1.gff','r')
    temp_gff = headle_gff.read()
    gene_gff.write(temp_gff)
    gene_gff.close()
    headle_gff.close()
    os.remove('gblast.out_1.1c_2.3_s1_0_16_1')
    os.remove(f'{single_seq}_{os.path.basename(genome_file)}.blast')
    os.remove(f'{single_seq}_{os.path.basename(genome_file)}.blast.report')
    os.remove('perform.txt')
    
def reconsitution_gff(file_gff):
    gene_gff = open(f'{file_gff}','r')
    clean_line = []
    for line in gene_gff:
        if line[0] != '#':
            clean_line.append(line.split())
            
            for line_list in clean_line:
                if line_list[2].strip() == 'transcript':
                    line_list[2] = 'gene'
                elif line_list[2].strip() == 'coding_exon':
                    line_list[2] = 'exon'                
    gene_gff.close()
    
    id_list = []
    for line_list in clean_line:
        if line_list[2].strip() == 'gene':
            id_list.append(line_list[8].split(';')[0][3:])
    id_list.reverse()
    cluster_list = []
    while len(id_list) > 0:
        id_name = id_list.pop()
        temp_list=[]
        for line_list in clean_line:
            if line_list[8].find(id_name) > 0:
                temp_list.append(line_list)
        cluster_list.append(temp_list)

    for cluster in cluster_list:
        for line_list in cluster:
            if line_list[2] == 'gene':
                id = str(f'{line_list[8].split(";")[0][3:]}')
                name = str(line_list[8].split(';')[1][5:])
                line_list[8] = f'ID={id};Name={name}'
            elif line_list[2] == 'exon':
                line_list[8] = f'ID={id}.exon{cluster.index(line_list)};Parent=mrna.{id}'
    cluster_list.sort()
    
    for cluster in cluster_list:
        for line_list in cluster:
            if line_list[2] == 'gene':
                id = str(f'{line_list[8].split(";")[0][3:]}')
                add_list = line_list.copy()
                add_list[2] = 'mRNA'
                add_list[8] = f'ID=mrna.{id};Parent={id}'
                add_list[5] = '.'
                cluster.insert(cluster.index(line_list)+1,add_list)
            elif line_list[2] == 'exon':
                add_list = line_list.copy()
                add_list[2] = 'CDS'
                add_list[8] = f'ID=cds.{id};Parent=mrna.{id}'
                cluster.insert(cluster.index(line_list)+1,add_list)
    return cluster_list
    
alignscore = open('alignscore.txt','a')
score = open(f"{sys.path[0]}/alignscore.txt",'r')
alignscore.writelines(score)
alignscore.close()
score.close()

id_list = get_single_seq(args.input) 
for id in id_list:
    run_gblast(f'{id}.fa',args.genome)
    os.remove(f'{id}.fa') 
    

gff_file = open(args.output,'a')
gff = csv.writer(gff_file,delimiter='\t')
for cluster in reconsitution_gff('./gene_gff'):
    for line_list in cluster:
        gff.writerow(line_list)
gff_file.close()
os.remove('./gene_gff')
os.remove('./gblast.out_1.1c_2.3_s1_0_16_1.gff')
os.remove('./alignscore.txt')

print('-------completed-------')

