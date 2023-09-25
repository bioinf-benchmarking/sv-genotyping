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
        "cp local_data/hprc.vcf.gz {output}"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"
        #"wget -O {output} https://zenodo.org/record/7669083/files/cactus_filtered_ids_chm13.vcf.gz?download=1"


#hprc needs no filtering
rule get_hprc_variants:
    input:
        RawHPRCVariants.path()
    output:
        RealVariantSource.path(database_name="hprc", variant_type="all", file_ending="/filtered.vcf.gz")
    shell:
        "cp {input} {output}"


