#!/usr/bin/env python3
# coding: utf-8

"""
2020/08/26
author:guisen chen
email:thecgs001@foxmail.com
description: By genblastG software, homologous gene of target 
             species was predicted based on homologous protein 
             sequences, and GFF files was generated.
usage:./run_genblastG.py seq.fa genome.fa
"""

import os
import sys
import csv
import subprocess
from Bio import SeqIO

os.chdir('/home/guisen/temp/genblast') # 设置工作路径
print(f'当前工作目录为：{os.getcwd()}')

def get_single_seq(seq_file):
    """
    将multiple_seq割成以seq_id.fa为文件名的single_seq文件,并返回一个以seq_id为元素的list
    """
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
    """
    运行gblastG,并生成gene_cds，gene_pep，gene_gff文件，此软件的当前目录必须有alignscore.txt与blastall两个文件
    """
    cline = subprocess.Popen(f'genblastG -P blast -q {single_seq} -t {genome_file} -e 1e-10 -r 5 -c 0.8 -gff -o gblast.out',env = None,shell=True)
    cline.wait()
    gene_gff = open('gene_gff','a')
    headle_gff = open('gblast.out_1.1c_2.3_s1_0_16_1.gff','r')
    temp_gff = headle_gff.read()
    gene_gff.write(temp_gff)
    gene_gff.close()
    headle_gff.close()
    os.remove('gblast.out_1.1c_2.3_s1_0_16_1')   # 此文件多序列文件运行会很大，所以每次运行一个序列后都删除
    os.remove('gblast.out_1.1c_2.3_s1_0_16_1.gff')
    os.remove(f'{single_seq}_{os.path.basename(genome_file)}.blast')
    os.remove(f'{single_seq}_{os.path.basename(genome_file)}.blast.report')
    os.remove('perform.txt')
    
def reconsitution_gff(file_gff):
    """
    清洗file.gff
    """
    gene_gff = open(f'{file_gff}','r')
    
    #删除注释行，将gff_file的第3列的transcript改为mRNA，coding_exon改为CDS，并返回一个以行为元素的clean_list
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
    
    #返回一个以gene_id为元素的id_list,用于cluster
    id_list = []
    for line_list in clean_line:
        if line_list[2].strip() == 'gene':
            id_list.append(line_list[8].split(';')[0][3:])
    id_list.reverse()
    
    #返回以一个三重嵌套list,第一层为gene_structure_cluster
    cluster_list = []
    while len(id_list) > 0:
        id_name = id_list.pop()
        temp_list=[]
        for line_list in clean_line:
            if line_list[8].find(id_name) > 0:
                temp_list.append(line_list)
        cluster_list.append(temp_list)

    #改为标准的gff3格式，并对cluster以染色体顺序排序,返回cluster_list
    for cluster in cluster_list:
        for line_list in cluster:
            #print(line_list)
            if line_list[2] == 'gene':
                id = str(f'{line_list[8].split(";")[0][3:]}')
                name = str(line_list[8].split(';')[1][5:])
                line_list[8] = f'ID={id};Name={name}'
            elif line_list[2] == 'exon':
                line_list[8] = f'ID={id}.exon{cluster.index(line_list)};Parent={id}'
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
                add_list[8] = f'ID=cds.{id};Parent={id}'
                cluster.insert(cluster.index(line_list)+1,add_list)
    return cluster_list
    
#执行
id_list = get_single_seq(sys.argv[1]) 
for id in id_list:
    run_gblast(f'{id}.fa',sys.argv[2])
    os.remove(f'{id}.fa') # 删除生成的single_seq_file
    
#写入gff文件
gff_file = open('./gene.gff','a')
gff = csv.writer(gff_file,delimiter='\t')
for cluster in reconsitution_gff('./gene_gff'):
    for line_list in cluster:
        gff.writerow(line_list)
gff_file.close()
os.remove('./gene_gff')
print('-------completed-------')

