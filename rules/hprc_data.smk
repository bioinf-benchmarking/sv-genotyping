"""
Rules for downloading and preprocessing HPRC data
"""
from typing import Literal
from snakehelp import parameters


@parameters
class RawHPRCVariants:
    folder: Literal["real_data/hprc"]
    file_ending = "/variants.vcf.gz"


rule downlad_raw_hprc_variants:
    output:
        protected(RawHPRCVariants.path())
    shell:
        "cp local_data/hprc_biallelic.vcf.gz {output}"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"


#hprc needs no filtering
rule get_hprc_variants:
    input:
        vcf = RawHPRCVariants.path(),
        index = RawHPRCVariants.path(file_ending="/variants.vcf.gz.tbi")
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


