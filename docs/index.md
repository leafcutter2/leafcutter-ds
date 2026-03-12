# leafcutter-ds

A Python re-implementation of the [leafcutter](https://github.com/davidaknowles/leafcutter) differential splicing algorithm.

![Leafcutter logo](https://github.com/leafcutter2/leafcutter-ds/images/logo.png){width=200}
## What does leafcutter do?

**Annotation-free quantification of RNA splicing.**

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include:

- Easy detection of novel introns
- Modeling of more complex splicing events than exonic PSI
- Avoiding the challenge of isoform abundance estimation
- Simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

Original method by [Yang I. Li](https://thelilab.com/)<sup>1</sup>, [David A. Knowles](https://daklab.github.io/)<sup>1</sup>, [Jack Humphrey](https://jackhump.github.io/), Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, [Jonathan K. Pritchard](http://web.stanford.edu/group/pritchardlab/home.html).

<sup>1</sup> *Equal contribution*

See the [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/044107) and [Nature Genetics publication](https://www.nature.com/articles/s41588-017-0004-9) for details.

Python re-implementation by Scott I. Adamson and David A. Knowles.

## Why another implementation?

The original leafcutter uses RStan, which has historically been difficult to install for many users. This Python/Pyro implementation overcomes that dependency. It is also compatible with the new [leafcutter2](https://github.com/leafcutter2/leafcutter2) implementation, which annotates unproductive splicing events based on leafcutter clusters (see the [Nature Genetics publication](https://doi.org/10.1038/s41588-024-01872-x)).

## Quick start

```bash
pip install leafcutter
```

Or install from source:

```bash
git clone https://github.com/leafcutter2/leafcutter-ds.git
cd leafcutter-ds
pip install -e .
```

The installation provides five command-line tools:

- `leafcutter-cluster` — cluster introns from junction files
- `leafcutter-ds` — differential splicing analysis
- `leafcutter-bayes` — Bayesian differential splicing
- `leafcutter-gtf-to-exons` — Exon table preparation from gtf file for use with leafcutter-ds
- `leafcutter-prepare-phenotype` — Splicing phenotype table for downstream sQTL calling tools.

## Community

If you have usage questions, join the [Google group](https://groups.google.com/forum/#!forum/leafcutter-users), or creating an issue on the [github](https://github.com/leafcutter2/leafcutter-ds/issues).
