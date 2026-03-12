# Getting Started

## Installation

Install from PyPI:

```bash
pip install leafcutter
```

Or install from source for development:

```bash
git clone https://github.com/leafcutter2/leafcutter-ds.git
cd leafcutter-ds
pip install -e .
```

## Dependencies

Key dependencies installed automatically through pip:

- **numpy**, **pandas**, **scipy** — core data manipulation
- **pyro-ppl** / **torch** — Bayesian inference backend (for `leafcutter-bayes`)

## CLI Tools

### `leafcutter-cluster`

Clusters introns from RNA-seq junction files.

```bash
leafcutter-cluster --help
```

### `leafcutter-ds`

Runs differential splicing analysis between two groups.

```bash
leafcutter-ds --help
```

### `leafcutter-bayes`

Bayesian differential splicing using Pyro.

```bash
leafcutter-bayes --help
```

### `leafcutter-gtf-to-exons`

A utility tool to make the exons file used in leafcutter-ds from a gtf file.

```bash
leafcutter-gtf-to-exons --help
```

### `leafcutter-prepare-phenotype`

A tool to prepare a splicing phenotype table for downstream sQTL calling tools.

```bash
leafcutter-prepare-phenotype --help
```

