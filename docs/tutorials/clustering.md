# Clustering Tutorial

This tutorial walks through intron clustering from RNA-seq junction files using `leafcutter-cluster`.

## Prerequisites

- Aligned RNA-seq BAM files (or pre-extracted junction files)
- A file listing junction file paths (one per line)

## Intron clustering options:

- Clustering with [leafcutter2](https://github.com/leafcutter2/leafcutter2), which annotates unproductive splicing events.
- Clustering with leafcutter-cluster available through this package outlined below.

## Step 1: Extract junctions

If starting from BAM files, extract splice junctions using [regtools](https://regtools.readthedocs.io/):

```bash
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 sample.bam -o sample.junc
```

Create a file listing all junction files:

```bash
ls *.junc > junc_files.txt
```

## Step 2: Run clustering

```bash
leafcutter-cluster \
    --juncfiles junc_files.txt \
    --outprefix my_study \
    --minclureads 30 \
    --mincluratio 0.001 \
    --maxintronlen 100000 \
    --rundir output/
```

<!-- TODO: add figure showing example cluster output -->

## Output description: leafcutter-cluster (this workflow)

### Counts file (`_perind_numers.counts.gz`)

The main counts file has one row per intron and one column per sample. The first column is the intron ID; subsequent columns are the counts corresponding to each intron junction count. This is the file passed to `leafcutter-ds`.

```
chrom                              sample1   sample2   ...
chr1:14829:14970:clu_1_-           348      351
chr1:14829:15021:clu_1_-           0        1
chr1:15947:16607:clu_2_-           13       3
```

| Field | Description |
|-------|-------------|
| `chr:start:end` | Intron coordinates (1-based, inclusive) |
| `clu_N` | Cluster ID |
| `strand` | `+`, `-`, or `NA` |
| `count` | Reads supporting this intron for that junction |

### Numerators file (`_perind.counts.gz`)

Same layout, but each sample contains the raw read count / the total read count the cluster for that sample. 

| File | Description |
|------|-------------|
| `my_study_perind_numers.counts.gz` | Per sample junction counts; input to differential splicing |
| `my_study_perind.counts.gz` | Per-sample intron usage ratios (numerator/total) |


## Options reference

| Option | Default | Description |
|--------|---------|-------------|
| `-j / --juncfiles` | *(required)* | Path to a text file listing junction files, one per line |
| `-o / --outprefix` | `leafcutter` | Prefix for all output file names |
| `-r / --rundir` | `./` | Directory to write output files |
| `-l / --maxintronlen` | `100000` | Maximum intron length in bp; longer junctions are discarded |
| `-m / --minclureads` | `30` | Minimum total reads across all samples for a cluster to be retained |
| `-p / --mincluratio` | `0.001` | Minimum fraction of cluster reads required to keep an individual intron |
| `-c / --cluster` | — | Provide a pre-computed refined cluster file to skip the clustering step |
| `-k / --nochromcheck` | `False` | Disable the check that chromosomes are named `chr1…chrY` or `1…Y` |
| `-C / --includeconst` | `False` | Also output constitutive introns (used in only one cluster) |
| `-q / --quiet` | `False` | Suppress progress messages |

## Output description: leafcutter2
See [leafcutter2](https://github.com/leafcutter2/leafcutter2) github for details. 
### Output:
- `leafcutter2.cluster_ratios.gz` a table quantifying splice junction read counts for each intron, divided by the total number of reads observed for the intron cluster to which the intron belongs. Each row corresponds to a splice junction, with the first column indicating the splice junction ID and all subsequent columns corresponding to each sample. The splice junction ID has the following format: `chr10:134786:179993:clu_1_+:PR`, indicating the chromosome, start and end of the splice junction, the LeafCutter intron cluster to which it belongs (in this case, `clu_1_+`), and a label indicating the splice junction's function: **PR** (productive/protein-coding), **UP** (unproductive), **NE** (ambiguous in their functional effect) or **IN** (intergenic).
- `leafcutter2.junction_counts.gz` a table similar to `leafcutter2.cluster_ratios.gz`, except it only keeps track of the splice junction read counts for each intron. This can be used downstream for `leafcutter-ds`.
- `clustering/` a directory containing files relevant for clustering and annotation.
    - `clustering/leafcutter2_clusters` contains the intron clusters. Useful for skipping the clustering step in repeated runs.
    - Other files documenting stats from the clustering and classification algorithm. Useful for debugging.

## Next steps

Proceed to [Differential Splicing](differential-splicing.md) to test for group differences.
