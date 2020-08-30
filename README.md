# **genblastG_extension : 一个基于genblastG软件，利用同源蛋白预测新基因的项目流程** 

## 1. get_ncbi_protien

**脚本：** get_ncbi_protien.py

**用法：** `./get_ncbi_protien.py id.txt`

**简介：** 输入一个含 NCBI protien id 的 id.txt 文件 ，返回一个含序列的 protien.fa 文件，并在屏幕输出查询到与未查询的序列 id 。

## 2. run_genblastG

**脚本：** run_formatdb_nucl.py

**用法：** ./run_formatdb_nucl.py genome.fa

**简介：** 利用 blast [v2.2.19] 的 formatdb 选项构建一个核酸数据库，构建的核酸数据库在输入基因组目录下，<br />

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 所以 Liunx 必须安装 blast 且版本号不超过 [v2.2.19] 。

<br />

**脚本：** /run_genblastG.py

**用法：** `./run_genblastG.py seq.fa genome.fa`

**简介：** 利用 genblastG [v1.39] 











<br />
<br />
<br />
<br /> 

## **genblastG_extension:** 
***a homology-based prediction gene process by genblastG***

**author:** *guisen chen*  |  **email:** *thecgs001@foxmail.com*

### get_single_seq.py

**Usage:** `./get_single_seq.py mulitple_seq.fa`

**Description:** To split a fasta file include of mulitple seq into a single_seq file, and return a list of seq_id.

**Note:** The Biopython package must be installed.

### run_formatdb_nucl.py

**Usage:** `./get_single_seq.py genome.fa`

**Description:** Creating a nucleic acid library 

**Note:** Based on the blast software.

### run_genblastG.py

**Usage:** `./run_genblastG.py seq.fa genome.fa`

**Description:** By genblastG software, homologous gene of target species was predicted based on homologous protein sequences, and GFF files was generated.

**Note:** Base on the genblast software, Biopython package.

### Base on [genblastG](http://genome.sfu.ca/genblast/download.html)
