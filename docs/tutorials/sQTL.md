# sQTLs

To perform splicing QTL (sQTL) analysis using LeafCutter you’ll first need to preprocess your RNA-seq data as described in Steps 1 and 2 under Differential Splicing. sQTL analysis is a little more involved than differential splicing analysis, but we provide a script scripts/prepare_phenotype_table.py intended to make this process a little easier. We assume you will use FastQTL for the sQTL mapping itself, but reformatting the output if you want to use another tool (e.g. MatrixEQTL ) should be reasonably straightforward. The script is pretty simple: a) calculate intron excision ratios b) filter out introns used in less than 40% of individuals or with almost no variation c) output these ratios as gzipped txt files along with a user-specified number of PCs.

Usage is 
`leafcutter-prepare-phenotype example/leafcutter_ds/Geuvadis_M_vs_F_perind.counts_sample.gz -p 10`

where `-p 10` specifies you want to calculate 10 PCs for FastQTL to use as covariates.

FastQTL needs tabix indices. To generate these you’ll need tabix and bgzip which you may have as part of samtools, if not they’re now part of htslib, see https://github.com/samtools/htslib for installation instructions. With these dependencies installed you can run the script created by and pointed to by the output of `leafcutter-prepare-phenotype`, e.g. `example/leafcutter_ds/Geuvadis_M_vs_F_perind.counts_sample.gz_prepare.sh`.

We assume you’ll run FastQTL separately for each chromosome: the files you’ll need will have names like `example/leafcutter_ds/Geuvadis_M_vs_F_perind_numers.counts_sample.gz.phen_chr1`. The PC file will be e.g. `example/leafcutter_ds/Geuvadis_M_vs_F_perind.counts_sample.gz.PCs`.

## Usage ##
`usage: prepare_phenotype_table.py [-h] [-p NPCS] [--chromosome_blacklist CHROMOSOME_BLACKLIST] counts_file`

| Option | Default | Description |
|--------|---------|-------------|
| counts_file | None| Intron usage counts file (input_perind.counts.gz), output from leafcutter-cluster |
| -h, --help | None| Shows help message|
| -p, --pcs | 50 | Number of PCs to output |
| --chromosome_blacklist | chrX and chrY | File of chromosomes to exclude from analysis, one per line. If not provided, defaults to blacklisting X and Y |

