from typing import Literal, Union
from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"
configfile: "config/plots.yaml"

workflow.use_conda = True

from snakehelp import parameters, result
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
    folder_name: Literal["real_variants"] = "real_variants"
    database_name: Literal["1000genomes", "hprc"] = "1000genomes"
    variant_type: Literal["snps_indels", "svs", "all"] = "snps_indels"
    file_ending = "/unfiltered_population.vcf.gz"


@parameters
class RealVariantSourceSingleChromosome:
    variants: RealVariantSource
    chromosome: str = "1"
    file_ending = "/variants.vcf.gz"


@parameters
class SimulatedPopulation:
    variants: SimulatedVariantSource
    allele_frequency: float = 0.3
    correlation: float = 0.8
    file_ending = "/unfiltered_population.vcf.gz"


@parameters
class RawPopulation:
    source: Union[RealVariantSource, SimulatedPopulation]
    file_ending = "/unfiltered_population.vcf.gz"


@parameters
class PopulationWithoutIndividual:
    population: RawPopulation
    individual_id: int = 1
    file_ending = "/population_without_individual.vcf.gz"


@parameters
class Individual:
    population: RawPopulation
    individual_id: int = 1
    file_ending = "/individual.vcf"


@parameters
class StratifiedIndividual:
    """
    An individual where variants are filtered on a stratification
    """
    individual: Individual
    stratification_variant_type: Literal["snps", "indels", "snps_indels", "svs", "all"] = "snps_indels"
    stratification_type: Literal["all", "easy", "low-mappability", "repeats", "other-difficult"] = "all"
    file_ending = "/individual.vcf"


@parameters
class FilteredIndividual:
    """Additional filter on a Stratified individual"""
    individual: StratifiedIndividual
    individual_filter: Literal["none", "only_variants_in_population"] = "none"
    file_ending = "/individual.vcf"


@parameters
class AlleleFrequencyFilteredPopulation:
    population: PopulationWithoutIndividual
    allele_frequency_svs: float = 0.1
    allele_frequency_snps_indels: float = 0.3
    file_ending = "/allele_frequency_filtered_population.vcf.gz"


@parameters
class FilteredPopulation:
    """Population filtered on allele frequency"""
    population: AlleleFrequencyFilteredPopulation
    #allele_frequency_svs: float = 0.1
    #allele_frequency_snps_indels: float = 0.3
    n_individuals: int = 10
    file_ending = "/filtered_population.vcf.gz"


@parameters
class Reads:
    # Reads are implicitly simulated from the individual that is removed from the population
    individual: Individual
    read_source: Literal["simulated", "real"] = "simulated"
    read_length: int = 150
    coverage: float = 10.0
    snp_error_rate: float = 0.001
    file_ending = "/reads.fq.gz"


@parameters
class RealRawReads:
    individual: Individual
    file_ending = "/real_reads.fq.gz"


@parameters
# Necessay to group reads and population so that GenotypeResults have one dependency
class ReadsAndFilteredPopulation:
    population: FilteredPopulation  # FilteredPopulation has reference to Individual
    read_source: Literal["simulated", "real"] = "simulated"
    read_length: int = 150
    coverage: float = 10.0
    snp_error_rate: float = 0.001
    file_ending = "/reads.fq.gz"


@parameters
class GenotypeResults:
    reads: ReadsAndFilteredPopulation
    method: Literal["pangenie", "kage", "kage_no_imputation", "kage_multiallelic", "kage_with_glimpse", "pangenie_multiallelic", "paragraph"] = "kage"
    n_threads: int = 4
    file_ending = "/genotypes.vcf"


@parameters
class StratifiedGenotypeResults:
    genotype_results: GenotypeResults
    stratification_variant_type: Literal["snps", "indels", "snps_indels", "svs", "all"] = "snps_indels"
    stratification_type: Literal["all", "easy", "low-mappability", "repeats", "other-difficult"] = "all"

    individual_filter: Literal["none", "only_variants_in_population"] = "none"
    file_ending = "/genotypes.vcf"


@parameters
class GenotypeDebug:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"
    file_ending = "/debug.txt"


@result
class GenotypeReport:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"

@result
class VcfEvalReport:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"

@result
class GenotypeRecall:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"

@result
class GenotypeOneMinusPrecision:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"

@result
class GenotypeF1Score:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"

@result
class GenotypeRuntime:
    genotype_results: StratifiedGenotypeResults
    limit_accuracy_to_variant_type: Literal["all", "snps", "indels", "snps_indels", "svs"] = "all"



@parameters
class GenomeStratification:
    genome_build: GenomeBuild
    stratification_type: Literal["all", "easy", "low-mappability", "repeats", "other-difficult"] = "all"
    file_ending = ".bed"



print(FilteredPopulation.path())


include: "rules/variant_simulation.smk"
include: "rules/population_simulation.smk"
include: "rules/read_simulation.smk"
include: "rules/pangenie.smk"
include: "rules/evaluation.smk"
include: "rules/kage.smk"
include: "rules/tests.smk"
include: "rules/thousand_genomes_data.smk"
include: "rules/hprc_data.smk"
include: "rules/bwa.smk"
include: "rules/paragraph.smk"
include: "rules/real_reads.smk"
include: "rules/glimpse.smk"

# for plotting
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")
#include: "/home/ivar/dev/mapping-benchmarking/rules/plotting.smk"



#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/reference_genome.smk"
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/read_simulation.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/mason.smk", branch="master")
#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")
#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/plotting.smk"

