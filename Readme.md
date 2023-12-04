
# Pipeline for benchmarking pangenomic structural variant genotyping

This is a Snakemake pipeline for benchmarking the performance of various approach for genptyping/calling structural variation using "pangenomes".

The aim of this pipeline is to make it easy to explore various approaches and how they affect accuracy, such as filtering of pangenomes, pangenome size, pangenome complexity, structural vairant type, etc.

The pipeline is designed so that experiments can be configurated using YAML specification, and small runs should be possible to run without a lot of computational resources.


## Overview of pipeline

The main structure is as follows:

Input to any analysis is usually two things:

1) Some sort of pangenome in vcf format is always the input. The pangenome should have variants and phased genotypes for a set of individuals
2) Some individual that we want to treat as a *truth*. We genotype that individual and compare the genotypes to this truth and then compute som accuacy measures. The pipeline supports randomly picking truth individuals from the input pangenome, or to use individuals from remote sources (e.g. GIAB).

The pipeline then works roughly as follows:

1) The input pangenome is optionally filtered using various techniques (see parameters below)
2) Reads are either simulated for the truth genome or optionally downloaded from a remote sources
3) One or more gentypers are run
4) Genotype calls are filtered or stratified in variuos ways
5) Accuracy measures are computed (currently Truvari is supported, as well as some more simple naive measure)

The aim is that the user can configure runs without thinking about how steps are connected, but rather by which parameters are relevant. Currently, these parameters are supported:

* `genome_build`: The reference genome our pangenome is based on, e.g. `hg38`. Will be automatically downloaded from UCSC.
* `size`: This parameter refers to a genomic region (specified in `config.yaml`) and lets you run experiments on smaller genomic regions, e.g a single chromosome. It only makes sense to run on smaller regions with simulated reads (not real reads)
* .... (to be documented)



## How to run






## How to contribute
TODO


