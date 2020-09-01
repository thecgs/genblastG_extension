#!/usr/bin/env python3
# coding: utf-8

"""
2020/09/01
author:guisen chen
email:thecgs001@foxmail.com
usage:./EVM_gff_sort.py file.gff
"""

import sys

def EVM_gff_sort(gff_file):
    gff_file = open(gff_file, 'r')
    gff_sort = open('EVM_sort.gff','a')
    gff = gff_file.read()
    
    gff_list = gff.strip().split('\n\n')
    
    gff_list.sort()
    for line in gff_list:
        gff_sort.write(f'{line}\n\n')
        
    gff_file.close()
    gff_sort.close()
        
#执行
EVM_gff_sort(sys.argv[1])
print('--------completed--------')
