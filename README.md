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

## Example:

```bash
$ cd ./example
$ python3 run_formatdb_nucl.py genome.fa
$ python3 run_genblastG.py -g genome.fa -i seq.fa -o out.gff3
$ head out.gff3
scaffold292     genBlastG       gene    151014  161317  77.4653 +       .       ID=NP_001092224.1-R2-2-A1;Name=NP_001092224.1
scaffold292     genBlastG       mRNA    151014  161317  .       +       .       ID=mrna.NP_001092224.1-R2-2-A1;Parent=NP_001092224.1-R2-2-A1
scaffold292     genBlastG       exon    151014  151347  .       +       .       ID=NP_001092224.1-R2-2-A1.exon1;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       CDS     151014  151347  .       +       .       ID=cds.NP_001092224.1-R2-2-A1;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       exon    159761  159975  .       +       .       ID=NP_001092224.1-R2-2-A1.exon2;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       CDS     159761  159975  .       +       .       ID=cds.NP_001092224.1-R2-2-A1;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       exon    160161  160280  .       +       .       ID=NP_001092224.1-R2-2-A1.exon3;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       CDS     160161  160280  .       +       .       ID=cds.NP_001092224.1-R2-2-A1;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       exon    161186  161317  .       +       .       ID=NP_001092224.1-R2-2-A1.exon4;Parent=mrna.NP_001092224.1-R2-2-A1
scaffold292     genBlastG       CDS     161186  161317  .       +       .       ID=cds.NP_001092224.1-R2-2-A1;Parent=mrna.NP_001092224.1-R2-2-A1
```

## **Note:**

The [biopython](https://biopython.org/) package must be installed.

The pipline is  base on  [genblastG](http://genome.sfu.ca/genblast/download.html) (v1.38) software and [genblastg_patch](https://github.com/epaule/genblastg_patch) patch.

The gff3 file obtained are redundant, redundancy filtering will be developed later.
