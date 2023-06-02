from typing import Literal
from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"
configfile: "config/plots.yaml"

workflow.use_conda = True

from snakehelp import parameters, result
from mapping_benchmarking.config import WholeGenomeReads, Individual
import snakehelp
snakehelp.set_data_folder("data/")

@parameters
class GenomeBuild:
    genome_build: Literal["sacCer3", "hg38"] = "sacCer3"


@parameters
class ReferenceGenome:
    genome_build: GenomeBuild
    file_ending = "/reference.fa"


@parameters
class BaseGenome:
    genome_build: GenomeBuild
    size: Literal["small", "medium", "large"] = "small"
    file_ending = "/reference.fa"


@parameters
class VariantSource:
    base_genome: BaseGenome
    snp_rate: float = 0.001
    small_indel_rate: float = 0.0001
    sv_indel_rate: float = 0.0001
    file_ending = "/variants.vcf"


@parameters
class Population:
    variants: VariantSource
    n_individuals: int = 10
    allele_frequency: float = 0.3
    correlation: float = 0.8
    file_ending = "/population.vcf"


@parameters
class PopulationWithoutIndividual:
    population: Population
    individual_id: str = "simulated_1"
    file_ending = "/population.vcf"


@parameters
class Individual:
    population: Population
    individual_id: str = "simulated_1"
    file_ending = "/individual.vcf"


@parameters
class Reads:
    # Reads are implicitly simulated from the individual that is removed from the population
    population_without_individual: PopulationWithoutIndividual
    read_length: int = 150
    coverage: float = 10.0
    file_ending = "/reads.fq.gz"


@parameters
class GenotypeResults:
    reads: Reads
    method: Literal["pangenie", "kage", "kage_no_imputation"] = "kage"
    n_threads: int = 1
    file_ending = "/genotypes.vcf"


@result
class GenotypeRecall:
    genotype_results: GenotypeResults
    variant_type: Literal["all", "snps", "small_indels", "svs"] = "all"

@result
class GenotypeOneMinusPrecision:
    genotype_results: GenotypeResults
    variant_type: Literal["all", "snps", "small_indels", "svs"] = "all"

@result
class GenotypeF1Score:
    genotype_results: GenotypeResults
    variant_type: Literal["all", "snps", "small_indels", "svs"] = "all"


print("Path", GenotypeRecall.path())

include: "rules/variant_simulation.smk"
include: "rules/population_simulation.smk"
include: "rules/read_simulation.smk"
include: "rules/pangenie.smk"
include: "rules/evaluation.smk"
include: "rules/kage.smk"
# for plotting
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")



#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/reference_genome.smk"
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/read_simulation.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/mason.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")
#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/plotting.smk"

