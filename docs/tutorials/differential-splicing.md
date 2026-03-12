# Differential Splicing Tutorial

This tutorial demonstrates differential splicing analysis with `leafcutter-ds` (frequentist) or `leafcutter-bayes` (Bayesian).

## Prerequisites

- Clustering output from `leafcutter-cluster` (see [Clustering Tutorial](clustering.md))
- A sample groups file
- Optional: exons file

## Prepare groups file

Create a tab-separated file assigning each sample to a group:

```
sample1  control
sample2  control
sample3  case
sample4  case
```

For continuous phenotypes (e.g. age, RIN), use numeric values in the second column. Additional columns specify covariates/confounders:

```
sample1  control  batch1  25
sample2  control  batch2  30
sample3  case     batch1  27
sample4  case     batch2  31
```

Numeric covariate columns are z-scored; string columns are one-hot encoded (first level dropped).

## Optional: exons file
This is a tsv file derived from a gtf file that has the following format:
| chr | start | end | strand | gene_name |
|-----|-------|-----|--------|-----------|
|chromosome|exon start|exon end| gene strand|gene name|

This can be prepared from a gtf using the following command:

`leafcutter-gtf-to-exons input_gtf output_file`

This supports gzip compressed input.

## Run differential splicing

```bash
leafcutter-ds \
    my_study_perind_numers.counts.gz \
    groups.txt \
    --output_prefix ds_results \
    --exon_file gencode_exons.txt.gz
```

## Interpret output


Both tools write two output files.

### `<prefix>_cluster_significance.txt`

One row per tested cluster:

| Column | Description |
|--------|-------------|
| `cluster` | Cluster ID in `chr:clu_N_strand` format |
| `status` | Test outcome: `tested`, `NotEnoughSamples`, `SingleJunc`, etc. |
| `loglr` | Log-likelihood ratio statistic |
| `df` | Degrees of freedom |
| `p` | Nominal p-value |
| `p.adjust` | BH-adjusted p-value across all clusters |
| `genes` | Gene name(s) overlapping this cluster (requires `--exon_file`) |
| `annotations` | Splicing consequence categories present in the cluster (leafcutter2 annotated input only) |

### `<prefix>_effect_sizes.txt`

One row per intron:

| Column | Description |
|--------|-------------|
| `intron` | Intron ID in `chr:start:end:clu_N_strand` format |
| `logef_<group>` | Log effect size for each non-baseline group |
| `psi_<baseline>` | Estimated percent-spliced-in (PSI) in the baseline group |
| `psi_<group>` | Estimated PSI in the comparison group |
| `deltapsi_<group>` | Difference in PSI (`psi_<group>` − `psi_<baseline>`) |

If using the leafcutter2 workflow, the intron will also have the predicted intron annotation appended to it (so `chr:start:end:clu_N_strand:annotation`), which can be any of the following:
- PR (productive/protein-coding)
- UP (unproductive)
- NE (ambiguous in their functional effect)
- IN (intergenic)


## Options reference

The following options apply to both `leafcutter-ds`.

### Positional arguments

| Argument | Description |
|----------|-------------|
| `counts_file` | Intron usage counts (`.txt` or `.txt.gz`), output of `leafcutter-cluster` |
| `groups_file` | Tab-separated: column 1 = sample name, column 2 = group/phenotype, optional further columns = covariates |

### Optional arguments

| Option | Default | Description |
|--------|---------|-------------|
| `-0 / --baseline_group` | `Control` | Reference group for categorical comparisons; other groups are tested against this |
| `-o / --output_prefix` | `leafcutter_ds` | Prefix for the two output files |
| `-e / --exon_file` | — | Exon annotation file (`chr`, `start`, `end`, `strand`, `gene_name`). Used only to label clusters with gene names |
| `-s / --max_cluster_size` | ∞ | Skip clusters with more introns than this |
| `-i / --min_samples_per_intron` | `5` | Discard introns with at least one read in fewer than this many samples |
| `-g / --min_samples_per_group` | `3` | For categorical comparisons: require this many samples per group to have ≥ `min_coverage` reads |
| `-c / --min_coverage` | `20` | Minimum total reads in a cluster for it to be tested |
| `-u / --min_unique_vals` | `10` | For continuous phenotypes: minimum number of distinct phenotype values after coverage filtering |
| `--init` | `brr` | Initialisation strategy: `brr` (Bayesian ridge regression), `rr` (ridge regression), `mult` (multinomial logistic), or `0` (zero) |
| `-p / --num_threads` | `1` | Number of parallel threads (`leafcutter-ds` only) |
| `--timeit` | `False` | Print a breakdown of time spent at each stage (for benchmarking) |

Proceed to [Visualization](visualization.md) to visualize splicing differences from `leafcutter-ds`.
