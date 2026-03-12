# leafcutter-ds #
## Installation and basic information ##
This is a python re-implementation of the leafcutter differential splicing algorithm (originally [here](https://github.com/davidaknowles/leafcutter)).

You should be able to install with `python -m pip install leafcutter`.

Alternatively, you can clone this repo and install it with pip locally:
 ```
 git clone https://github.com/leafcutter2/leafcutter-ds.git
 cd leafcutter-ds
 pip install -e .
 ```

The installation process will install `leafcutter-cluster` and `leafcutter-ds` as command line tools.

This is compatible with the new [leafcutter2 implmetation](https://github.com/leafcutter2/leafcutter2) which annotates unproductive splicing events based on leafcutter clusters. See the [Nature Genetics publication](https://doi.org/10.1038/s41588-024-01872-x) for details.

### What does leafcutter do? ####
<img src="https://github.com/davidaknowles/leafcutter/blob/master/docs/logo.png" width="200"> **Annotation-free quantification of RNA splicing.**

[Yang I. Li](https://thelilab.com/)<sup>1</sup>, [David A. Knowles](https://daklab.github.io/)<sup>1</sup>, [Jack Humphrey](https://jackhump.github.io/), Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, [Jonathan K. Pritchard](http://web.stanford.edu/group/pritchardlab/home.html)

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. The advantages of this approach include:
* easy detection of novel introns
* modeling of more complex splicing events than exonic PSI
* avoiding the challenge of isoform abundance estimation
* simple, computationally efficient algorithms scaling to 100s or even 1000s of samples

For details please see our [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/044107) and corresponding [Nature Genetics publication](https://www.nature.com/articles/s41588-017-0004-9).

If you have usage questions we've setup a Google group here: <https://groups.google.com/forum/#!forum/leafcutter-users>

We've developed a leafcutter [shiny](https://shiny.rstudio.com/) app for visualizing leafcutter results: you can view an example [here](https://leafcutter.shinyapps.io/leafviz/). This shows leafcutter differential splicing results for a comparison of 10 brain vs. 10 heart samples (5 male, 5 female in each group) from [GTEx](https://www.gtexportal.org/home/).

### Why another implementation? ###
The short answer is that the Rstan dependecies that were difficult to install for many users, so we made this python/pyro version to overcome that. At the same time, we took the opportunity to make it compatible with new leafcutter methods to annotate unproductive splicing events. 

A much more detailed write-up on this version please see our readthedocs.


