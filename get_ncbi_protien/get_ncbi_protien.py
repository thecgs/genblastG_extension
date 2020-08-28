#!/usr/bin/env python3
#coding: utf-8

"""
2020/08/28
author:guisen chen
email:thecgs001@foxmail.com
description:According to the ID number, the protein sequence 
            was queried in batches from the NCBI protein database.
usage:./get_ncbi_protien seq_id.txt
"""

import sys
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'thecgs001@foxmail.com'

def get_ncbi_protien(seq_id):
    infile = open(seq_id,'r')
    outfile = open('protein.fa','a')
    error_list = []
    correct_list = []
    for seq_id in infile:
        try:
            return_query = Entrez.efetch(db='protein',rettype='gb',retmode='text',id=seq_id.strip())
            seq_record = SeqIO.read(return_query,'gb')
            outfile.write(f'>{seq_record.name}|{seq_record.description}\n{seq_record.seq}\n')
            correct_list.append(seq_id.strip())
        except:
            error_list.append(seq_id.strip())
    infile.close()
    outfile.close()
    return print(f'No result: {error_list}\nResult: {correct_list}')#返回查询和未查询到的id号

#执行    
get_ncbi_protien(sys.argv[1])

