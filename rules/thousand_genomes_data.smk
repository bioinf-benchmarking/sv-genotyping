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
        RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/raw.vcf.gz")
    output:
        RealVariantSource.path(database_name="1000genomes", variant_type="svs", file_ending="/raw_without_unknown_sequences.vcf.gz")
    conda:
        "../envs/filter_svs.yml"
    shell:
        """
        python scripts/filter_svs.py {input} | bgzip -c > {output}
        """


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


def all_files(wildcarsd):
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
        RealVariantSource.path(database_name="1000genomes", variant_type="all"),
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
    output:
        vcf = AlleleFrequencyFilteredPopulation.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        "python scripts/filter_vcf_on_allele_frequency.py {input.vcf} {wildcards.allele_frequency_svs} {wildcards.allele_frequency_snps_indels} | bgzip -c > {output.vcf}"


rule filter_population:
    """For filtering any population (real or simulated) on number of individuals and allele frequency"""
    input:
        vcf = AlleleFrequencyFilteredPopulation.path(),
        index = AlleleFrequencyFilteredPopulation.path(file_ending="/allele_frequency_filtered_population.vcf.gz.tbi"),
        sample_names = PopulationWithoutIndividual.path(file_ending="/population_without_individual.{n_individuals}_random_sample_names.txt")
    output:
        vcf = FilteredPopulation.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view -S {input.sample_names} -O z -o {output.vcf} {input.vcf}"



