"""
Rules for downloading 1000 genomes SNPs, indels and SVs and combining them into one dataset
"""

from snakehelp import parameters
from typing import Literal

@parameters
class Raw1000GenomesVariants:
    folder: Literal["real_data/1000genomes_snps_indels"]
    chromosome: str = "2"
    file_ending = "/variants.vcf.gz"

@parameters
class Raw1000GenomesStructuralVariants:
    folder: Literal["real_data/1000genomes_svs"]
    file_ending = "/variants.vcf.gz"

@parameters
class Filtered1000GenomesStructuralVariants:
    variants: Raw1000GenomesStructuralVariants
    file_ending = "/filtered/variants.vcf.gz"

@parameters
class Filtered1000GenomesVariants:
    variants: Raw1000GenomesVariants
    allele_frequency: float = 0.01
    file_ending = "/variants.vcf.gz"


@parameters
class FilteredOnSamples1000GenomesVariants:
    # filtered on individuals in the SV dataset
    variants: Filtered1000GenomesVariants
    file_ending = "/sample_filtered/variants.vcf.gz"


@parameters
class FilteredOnSamples1000GenomesStructuralVariants:
    variants: Filtered1000GenomesStructuralVariants
    file_ending = "/filtered/sample_filtered/variants.vcf.gz"

@parameters
class FilteredMerged1000GenomesVariants:
    folder: Literal["real_data/1000genomes"]
    allele_frequency: float = 0.01
    file_ending = "/variants.vcf.gz"

print(FilteredMerged1000GenomesVariants.path())


rule download_thousand_genomes_snps_indels_for_chromosome:
    output:
        RealVariantSourceSingleChromosome.path(database_name="1000genomes", variant_type="snps_indels")
    params:
        #url = lambda wildcards: f"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{wildcards.chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
        url = lambda wildcards: f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{wildcards.chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    shell:
        # tmp for testing, keep only few lines
        #"true || curl -NL {params.url} 2>/dev/null  "
        "curl -NL {params.url}"
        "| zcat "
        #"| head -n 10000 "
        "| python scripts/filter_variants_with_n.py "
        #" > {output}.tmp && "
        #" && python scripts/remove_overlapping_indels.py {output}.tmp  "
        "| bgzip -c > {output}  "
        " || true"


rule filter_variants_on_allele_frequency:
    input:
        Raw1000GenomesVariants.path()
    output:
        Filtered1000GenomesVariants.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view -q {wildcards.allele_frequency} {input} -Oz -o {output}"


rule download_1000genomes_svs:
    output:
        RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/raw.vcf.gz")
    shell:
        "wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz"



#only keeping SVs with known sequences
rule filter_1000genomes_svs:
    input:
        variants = RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/raw.vcf.gz"),
        reference = ReferenceGenome.path()
    output:
        RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/raw_without_unknown_sequences.vcf.gz")
    #conda:
    #    "../envs/bcftools.yml"
    shell:
        #"python scripts/filter_svs.py {input} | bgzip -c > {output}"
        "kage preprocess_sv_vcf -v {input.variants} -f {input.reference} | bgzip -c > {output}"



rule subset_variants_on_dataset_chromosomes:
    input:
        vcf = RealVariantSource.path(variant_type="svs", file_ending="/raw_without_unknown_sequences.vcf.gz"),
        index = RealVariantSource.path(variant_type="svs", file_ending="/raw_without_unknown_sequences.vcf.gz.tbi"),
    output:
        RealVariantSource.path(variant_type="svs", file_ending="/unfiltered_population.vcf.gz")  # may need variant_type=svs to not conflict with snps/indels
    params:
        chromosomes = lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"]
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --regions {params.chromosomes} {input.vcf} -Oz -o {output}
        """




rule subset_on_individuals_snps_indels:
    input:
        variants=Filtered1000GenomesVariants.path(),
        individiuals="data/real_data/1000genomes_snps_indels/samples_common.txt"
    output:
        FilteredOnSamples1000GenomesVariants.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
            bcftools view --samples-file {input.individiuals} {input.variants} -Oz -o {output}
            """



def snps_indels_files(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"].split(",")
    return [RealVariantSourceSingleChromosome.path(database_name="1000genomes", variant_type="snps_indels", chromosome=chromosome) for chromosome in chromosomes]


def all_files(wildcards):
    chromosomes = config["real_datasets"]["1000genomes"]["chromosomes"].split(",")
    print(chromosomes)
    return [FilteredOnSamples1000GenomesVariants.path(chromosome=chromosome) for chromosome in chromosomes] + \
            [FilteredOnSamples1000GenomesStructuralVariants.path()]

# not used??

rule merge_1000genomes_snps_indels:
    input:
        snps_indels_files
    output:
        RealVariantSource.path(database_name="1000genomes", variant_type="snps_indels", file_ending="/unfiltered_population.vcf.gz")
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools concat --naive {input} -Oz -o {output}"


# snps_indels is only for snps indels, these vcfs needs to be merged
# svs and "all" are generally  from one single vcf which needs to be subset
ruleorder: merge_1000genomes_snps_indels > subset_variants_on_dataset_chromosomes

rule tabix:
    input:
        "{variants}.vcf.gz"
    output:
        "{variants}.vcf.gz.tbi"
    conda:
        "../envs/bgzip.yml"
    shell:
        "tabix -p vcf {input}"


rule merge_1000genomes_variant_types:
    input:
        snps_indels = RealVariantSource.path(database_name="1000genomes", variant_type="snps_indels", file_ending="/unfiltered_population.vcf.gz"),
        snps_indels_index = RealVariantSource.path(database_name="1000genomes", variant_type="snps_indels", file_ending="/unfiltered_population.vcf.gz.tbi"),
        svs = RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/unfiltered_population.vcf.gz"),
        svs_index = RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/unfiltered_population.vcf.gz.tbi"),
    output:
        RealVariantSource.path(database_name="1000genomes", variant_type="all", file_ending="/unfiltered_population_all.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        # first find individuals that are shared between the two
        """
        bcftools query -l {input.snps_indels} > {input.snps_indels}.samples &&
        bcftools query -l {input.svs} > {input.svs}.samples &&
        comm -12 {input.snps_indels}.samples {input.svs}.samples > {output}.samples &&
        bcftools view --samples-file {output}.samples {input.snps_indels} -Oz -o {input.snps_indels}.subset.vcf.gz &&
        bcftools view --samples-file {output}.samples {input.svs} -Oz -o {input.svs}.subset.vcf.gz &&
        tabix -p vcf {input.snps_indels}.subset.vcf.gz &&
        tabix -p vcf {input.svs}.subset.vcf.gz &&
        bcftools concat --allow-overlaps {input.snps_indels}.subset.vcf.gz {input.svs}.subset.vcf.gz -Oz -o {output}
        """

rule remove_snps_indels_inside_svs:
    input:
        RealVariantSource.path(database_name="1000genomes",variant_type="all",file_ending="/unfiltered_population_all.vcf.gz"),
    output:
        RealVariantSource.path(database_name="1000genomes",variant_type="all",file_ending="/unfiltered_population.vcf.gz")
    shell:
        """
        kage filter_snps_indels_covered_by_svs --vcf {input} -l 50 | bgzip -c > {output}
        """


#hacky rule to filter out a dataset with only snps and indels and give it a database name
#rule get_1000genomes_only_snps_and_indels_dataset:
#    input:
#        snps_indels = RealVariantSource.path(database_name="1000genomes", variant_type="snps_indels"),
#    output:
#        snps_indels = RealVariantSource.path(database_name="1000genomes_only_snps_indels", variant_type="snps_indels"),
#    shell:
#        "cp {input.snps_indels} {output.snps_indels}"


rule get_all_sample_names_from_vcf:
    input:
        vcf="{variants}.vcf.gz",
        index="{variants}.vcf.gz.tbi"
    output:
        "{variants}.all_sample_names.txt"
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools query -l {input.vcf} > {output}"


rule get_all_sample_names_from_vcf_nogz:
    input:
        "{variants}.vcf"
    output:
        "{variants}.all_sample_names.txt"
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools query -l {input} > {output}"


rule get_random_sample_names_from_vcf:
    input:
        "{variants}.all_sample_names.txt"
    output:
        out="{variants}.{n}_random_sample_names.txt"
    conda:
        "../envs/python.yml"
    shell:
        "python scripts/shuffle_lines.py {input} {config[random_seed]} {wildcards.n} {output}"


rule get_random_sample_name_from_number:
    input:
        "{variants}.10000_random_sample_names.txt"  # 10000 to get all sample names as input
    output:
        "{variants}.random_sample_number_{n}.txt"
    shell:
        "head -n {wildcards.n} {input} | tail -n 1 > {output}"


rule filter_population_allele_frequency:
    input:
        vcf= PopulationWithoutIndividual.path(),
        index=PopulationWithoutIndividual.path(file_ending="/population_without_individual.vcf.gz.tbi"),
        ref = ReferenceGenome.path()
    output:
        vcf = AlleleFrequencyFilteredPopulation.path()
    #conda:
    #    "../envs/bcftools.yml"
    shell:
        "python scripts/filter_vcf_on_allele_frequency.py {input.vcf} {wildcards.allele_frequency_svs} {wildcards.allele_frequency_snps_indels} | bgzip -c > {output.vcf}.tmp.vcf.gz && "
        "kage filter_low_freq_alleles -r {input.ref} -f {wildcards.minor_allele_in_multiallelic_variants_min_frequency} -v {output.vcf}.tmp.vcf.gz --only-deletions True | bgzip -c > {output.vcf} "


rule filter_population:
    """For filtering any population (real or simulated) on number of individuals and allele frequency"""
    input:
        vcf = AlleleFrequencyFilteredPopulation.path(),
        index = AlleleFrequencyFilteredPopulation.path(file_ending="/allele_frequency_filtered_population.vcf.gz.tbi"),
        sample_names = PopulationWithoutIndividual.path(file_ending="/population_without_individual.{n_individuals}_random_sample_names.txt")
    output:
        vcf = FilteredPopulation.path()
    params:
        # if we only want variants supported by individuals, filter so that min allele count is 1
        filter_command = lambda wildcards: " | bcftools view --min-ac 1  " if wildcards.population_type == "only_variants_supported_by_individuals" else ""
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view -S {input.sample_names}  {input.vcf} {params.filter_command} | bgzip -c > {output}"


rule collapse_similar_alleles_in_population:
    input:
        vcf = FilteredPopulation.path(),
        index = FilteredPopulation.path(file_ending="/filtered_population.vcf.gz.tbi"),
    output:
        vcf = FilteredPopulationWithCollapsedAlleles.path()
    conda:
        "../envs/truvari.yml"
    shell:
        "truvari collapse -r 1 -p {wildcards.collapse_threshold} -P {wildcards.collapse_threshold} -i {input.vcf} | bgzip -c > {output}"


rule collapse_similar_alleles_in_population_no_collapse:
    # When similarity parameter is 1.0, no collapsing
    input:
        vcf=FilteredPopulation.path(),
        index=FilteredPopulation.path(file_ending="/filtered_population.vcf.gz.tbi"),
    output:
        vcf = FilteredPopulationWithCollapsedAlleles.path(collapse_threshold=1.0)
    conda:
        "../envs/truvari.yml"
    shell:
        "cp {input.vcf} {output.vcf}"


ruleorder: collapse_similar_alleles_in_population_no_collapse > collapse_similar_alleles_in_population


