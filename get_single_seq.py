#!/usr/bin/env python3
# coding: utf-8

import sys
import os
from Bio import SeqIO

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
print(get_single_seq(sys.argv[1]))

# 陈桂森 2020/08/25

