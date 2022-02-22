# **genblastG_extension**

**A pipline of homology-based prediction gene**

**Author:** guisen chen  |  **Email:** <thecgs001@foxmail.com>

## Step 1:  Installation

```bash
$ git clone git@github.com:thecgs/genblastG_extension.git
$ export PATH=$PATH:user/genblastG_extension >> ~/.bashrc
$ source ~/.bashrc
```

## Step 2: creating a nucleic acid database

```bash
$ python3 run_formatdb_nucl.py genome.fa
```

## Step 3:  homology-based prediction

```bash
$ python3 run_genblastG.py -g genome.fa -i seq.fa -o out.gff3
```

## **Note:**

The [biopython](https://biopython.org/) package must be installed.

The pipline is  base on  [genblastG](http://genome.sfu.ca/genblast/download.html) (v1.38) software and [genblastg_patch](https://github.com/epaule/genblastg_patch) patch.

The gff3 file obtained are redundant, Redundancy filtering will be developed later.
