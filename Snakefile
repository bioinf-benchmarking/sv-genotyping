from typing import Literal, Union
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
    genome_build: Literal["sacCer3", "hg38", "chm13"] = "sacCer3"


@parameters
class ReferenceGenome:
    genome_build: GenomeBuild
    file_ending = "/reference.fa"


@parameters
class BaseGenome:
    genome_build: GenomeBuild
    size: str = "small"
    file_ending = "/reference.fa"


@parameters
class SimulatedVariantSource:
    base_genome: BaseGenome
    snp_rate: float = 0.001
    small_indel_rate: float = 0.0001
    sv_indel_rate: float = 0.0001
    file_ending = "/variants.vcf"


@parameters
class RealVariantSource:
    base_genome: BaseGenome
    folder_name: Literal["real_variants"]
    database_name: Literal["1000genomes", "hprc"] = "1000genomes"
    variant_type: Literal["snps_indels", "svs", "all"] = "snps_indels"
    file_ending = "/variants.vcf.gz"


@parameters
class RealVariantSourceSingleChromosome:
    variants: RealVariantSource
    chromosome: str = "1"
    file_ending = "/variants.vcf.gz"


@parameters
class SimulatedPopulation:
    variants: SimulatedVariantSource
    n_individuals: int = 10
    allele_frequency: float = 0.3
    correlation: float = 0.8
    file_ending = "/population.vcf"



@parameters
class FilteredRealPopulation:
    population: RealVariantSource
    allele_frequency_svs: float = 0.1
    allele_frequency_snps_indels: float = 0.3
    n_individuals: int = 10
    file_ending = "/population.vcf"

@parameters
class Population:
    population: Union[FilteredRealPopulation, SimulatedPopulation]
    file_ending = "/population.vcf"


@parameters
class PopulationWithoutIndividual:
    population: Population
    individual_id: int = 1
    file_ending = "/population.vcf"


@parameters
class Individual:
    population: Population
    individual_id: int = 1
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
    method: Literal["pangenie", "kage", "kage_no_imputation", "kage_multiallelic", "pangenie_multiallelic"] = "kage"
    n_threads: int = 1
    file_ending = "/genotypes.vcf"


@parameters
class GenotypeDebug:
    genotype_results: GenotypeResults
    variant_type: Literal["all", "snps", "small_indels", "svs"] = "all"
    file_ending = "/debug.txt"


@result
class GenotypeReport:
    genotype_results: GenotypeResults
    variant_type: Literal["all", "snps", "small_indels", "svs"] = "all"

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


print("Path", Population.path())

include: "rules/variant_simulation.smk"
include: "rules/population_simulation.smk"
include: "rules/read_simulation.smk"
include: "rules/pangenie.smk"
include: "rules/evaluation.smk"
include: "rules/kage.smk"
include: "rules/tests.smk"
include: "rules/thousand_genomes_data.smk"
include: "rules/hprc_data.smk"
# for plotting
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")



#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/reference_genome.smk"
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/read_simulation.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/mason.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")
#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/plotting.smk"

