## Python-script

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
