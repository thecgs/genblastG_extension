#!/usr/bin/env python3
#coding: utf-8

"""
2020/08/26
author:guisen chen
email:thecgs001@foxmail.com
description:creating a nucleic acid library
usage:./run_formatdb_nucl.py genome.fa
"""

import os
import sys
import subprocess

def run_formatdb_nucl(genome_file):
    """
    运行blast程序的formatdb，创建一个核酸库
    """
    cline = subprocess.Popen(f"{sys.path[0]}/formatdb -i {genome_file} -p F -o T",env=None,shell=True)
    cline.wait()
    cline = subprocess.Popen(f"mv {genome_file}.* ./",env=None,shell=True)
    cline.wait()
    cline = subprocess.Popen(f"cp {genome_file} ./",env=None,shell=True)
    cline.wait()
    
    handle = open(f'formatdb.log','r')
    formatdb_log = handle.read()
    handle.close()
    os.remove(f'formatdb.log')
    return print(formatdb_log)

run_formatdb_nucl(sys.argv[1])