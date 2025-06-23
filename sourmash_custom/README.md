# Sourmash Custom  
The sourmash_custom.nf subworkflow makes a custom sourmash database from a collection of fasta assembly files. 

## Expected Input
The sourmash_custom.nf workflow takes a text file with a single sequence per line, `db_manifest` (default=""). Each entry should contain the absolute path to the sequence in fasta format. Entries can be gzip compressed or a standard fasta input as per the requirements of sourmash.

Example manifest:
```
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/003/725/GCA_003003725.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/003/615/GCA_003003615.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/435/GCA_003524435.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/565/GCA_003524565.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/745/GCA_003524745.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/145/GCA_003524145.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/305/GCA_003524305.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/095/GCA_003524095.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/365/GCA_003524365.1_genomic.fna.gz
/data/pam/software/GTDB/release220/genomic_files_reps/gtdb_genomes_reps_r220/database/GCA/003/524/845/GCA_003524845.1_genomic.fna.gz
```

The `sourmash_database.config` file takes the kmer size, `klen` (default=31), and a scaling parameter, `sketch_size`  (default=10000). These are used in the sourmash parameter string e.g. `-p scaled=10000,k=31` to express the sketch sizes of the sequences within the database. Additional parameters include the database/index name, `db_name` (default="custom_db") adds a prefix to the final output file e.g. `custom_db_s10000k31_sourmash.zip`. 

## Expected Output
The pipeline outputs a `custom_db_s10000k31_sourmash.zip` file that contains all of the sketchs in the chosen output directory (default="./results/sourmash_db"). The contents of the database can be viewed using the `sourmash sig summarize` command.

```
sourmash sig summarize results/sourmash_db/custom_db_s10000k31_sourmash.zip 
```

The `sourmash sig summarize` command gives the version of sourmash used to generate the database as well as details about the the number of sketches, hashes, and sketch parameters. 

```
== This is sourmash version 4.5.0. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

** loading from 'results/sourmash_db/custom_db_s10000k31_sourmash.zip'
path filetype: ZipFileLinearIndex
location: /results/sourmash_db/custom_db_s10000k31_sourmash.zip
is database? yes
has manifest? yes
num signatures: 10
** examining manifest...
total hashes: 3213
summary of sketches:
   10 sketches with DNA, k=31, scaled=10000           3213 total hashes
```

In addition to the sourmash sketches, the pipeline uses the `parse_stb.py` script from dRep to create an scaffold-to-bin file for use in inStrain (default="custom_db_drep.stb").

Example scaffold-to-bin file:
'''
PVWD01000001.1	GCA_003003615.1_genomic.fna.gz
PVWD01000359.1	GCA_003003615.1_genomic.fna.gz
PVWD01000087.1	GCA_003003615.1_genomic.fna.gz
PVWD01000360.1	GCA_003003615.1_genomic.fna.gz
PVWD01000361.1	GCA_003003615.1_genomic.fna.gz
PVWD01000362.1	GCA_003003615.1_genomic.fna.gz
PVWD01000363.1	GCA_003003615.1_genomic.fna.gz
PVWD01000364.1	GCA_003003615.1_genomic.fna.gz
PVWD01000088.1	GCA_003003615.1_genomic.fna.gz
'''