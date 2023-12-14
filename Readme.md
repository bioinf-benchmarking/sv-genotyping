
# Pipeline for benchmarking pangenomic structural variant genotyping

This is a Snakemake pipeline for benchmarking the performance of various approaches for genotyping/calling structural variation using "pangenomes".

The aim of this pipeline is to make it easy to explore various methods and parameters and how they affect accuracy, e.g. how as filtering of pangenomes, pangenome size, pangenome complexity, structural vairant type, etc. affect accuracy.

The pipeline is designed so that experiments can be configurated using YAML specification, and small experiments should be possible to run without a lot of computational resources.


## Overview of pipeline

Input to any analysis is usually two things:

1) Some sort of pangenome in VCF format is always the input. The pangenome should have variants and phased genotypes for a set of individuals.
2) Some individual that we want to treat as a *truth*. We genotype that individual and compare the genotypes to this truth and then compute som accuacy measures. The pipeline supports randomly picking truth individuals from the input pangenome, or to use individuals from remote sources (e.g. GIAB).

The pipeline then works roughly as follows:

1) The input pangenome is optionally filtered using various techniques (see parameters below).
2) Reads are either simulated from the truth genome or optionally downloaded from a remote sources
3) One or more gentypers are run
4) Genotype calls are filtered or stratified in variuos ways
5) Accuracy measures are computed (currently Truvari is supported, as well as some more simple naive measure)
6) Some sort of plot is generated to visualize the results. This is done automatically based on YAML configuration.

The aim is that the user can configure runs without thinking about how steps are connected or how each step is performed, but rather by which parameters are relevant. Currently, these parameters are supported:

* `genome_build`: The reference genome our pangenome is based on, e.g. `hg38`. Will be automatically downloaded from UCSC.
* `size`: This parameter refers to a genomic region (specified in `config.yaml`) and lets you run experiments on smaller genomic regions, e.g a single chromosome. It only makes sense to run on smaller regions with simulated reads (not real reads)
* .... (to be documented)


## How to setup the pipeline and run small experiments locally
This assumes you have conda and Snakemake already installed.

### Step 1: Clone this repository
```bash
git clone 
```

### Step 2: Install Python requirements
```bash
cd sv-genotyping
pip install -r requirements.txt
```

### Step 3: Test that the pipeline works
```
snakemake test
```
If the following command does not give any error, the basic functionality should work, and you can proceed by running full experiments (see below).

## Preconfigured experiments

The pipeline contains variuous pre-configured experiments (defined in config/plots.yaml). Running these is straight forward given that you have enough computational resources:




## Configuring a custom experiment


### Creating a plot

This pipeline follows the Snakemake principles, meaning that the user defines what the final result should be, and then the pipeline tries to run the necessary jobs for creating that output. For instance, you can ask for a plot where the x-axis is something, the y-axis is something and the pipeline will try to run what is needed to generate that plot. "Something" needs to be a valid *parameter* or *result_type* (see above), and based on that, the pipeline figures out what rules to run. For instance, the x-axis could be `method` (i.e. genotyper) and the y-axis can be `f1-score` and the pipeline will then run all methods and capture the F1 score and present that.

You can add a plot type by adding a configuration under `plot_types` in `config/plots.yaml`. Example:

```yaml
plot_types:
  accuracy_vs_read_length:
    type: line
    x: read_length
    y: TruvariF1Score
    color: method
    facet_col: stratification_variant_type 
```

The above defines a plot type with the name `accuracy_vs_read_length`. We tell the pipeline to make a plot type where the x-axis is `read_length` and the y-axis is `TruvariF1Score` (which is a result type that is based on F1 score computed by Truvari). The "color" is `method`, which means that we want one line for each available method (color is the term that **Plotly** uses). We want to repeat this plot for different "variant types" along the columns (specified by `facet_col`). Note that only `x` and `y` are mandatory, the rest can be ommited (in that case only a single plot is created).

The above **only** specifies a **plot_type**, not a plot. For instance, we say that the x-axis should have `method`, but we don't tell it *which* methods to include. Using this plot type, we can specify actual plots with data, under the `plots`-section in the same YAML-file:

```yaml

plots:
    my_plot:
    plot_type: accuracy_vs_coverage
    parameters:
      coverage: [2, 10, 30]
      method: [kage_with_glimpse]
      sv_indel_rate: 0.01
    type_limits:
      source: SimulatedPopulation
  
```

All parameters thare are not specified will be set to their default values (specified in the parameter configuration in Snakefile). Since source is an ambiguous parameter (the source can be either simulated or a real variant source), we need to specify under `type_limits` what the source is.

We can now ask Snakemake to generate my_plot:

```bash
snakemake plots/my_plot.png
```
.. which should create something like this:

![Plot example](...)



