"""
Rules for downloading and preprocessing HPRC data
"""
from typing import Literal
from snakehelp import parameters


@parameters
class RawHPRCVariants:
    folder: Literal["real_data/hprc"]
    file_ending = "/variants.vcf.gz"


"""
Tmp rules to use local data instead of download.
Should be changed to download from zenodo
"""
rule downlad_raw_hprc_variants:
    output:
        protected(RawHPRCVariants.path(file_ending="/hprc_biallelic_chm13.vcf.gz"))
    shell:
        "cp local_data/hprc_biallelic_chm13.vcf.gz {output}"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"


rule downlad_raw_hprc_variants_hg38:
    output:
        protected(RawHPRCVariants.path(file_ending="/hprc_biallelic_hg38.vcf.gz"))
    conda:
        "../envs/bcftools.yml"
    shell:
        "zcat local_data/cactus_filtered_ids_biallelic.vcf.gz "
        " | sed 's/chr//g' | bgzip -c > {output}  "
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"


#hprc needs no filtering
rule get_hprc_variants:
    input:
        vcf = RawHPRCVariants.path(file_ending="/hprc_biallelic_{genome_build}.vcf.gz"),
        index = RawHPRCVariants.path(file_ending="/hprc_biallelic_{genome_build}.vcf.gz.tbi")
    output:
        RealVariantSource.path(database_name="hprc", variant_type="all", file_ending="/unfiltered_population.vcf.gz")
    params:
        chromosomes = lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"]
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --regions {params.chromosomes} {input.vcf} -Oz -o {output}
        """


