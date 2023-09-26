


rule samtools_index:
    input:
        "{sample}.fa",
    output:
        "{sample}.fa.fai",
    wrapper:
        "v1.21.2/bio/samtools/faidx"
    
    
def get_reference_url(wildcards):
    if wildcards.genome_build == "chm13":
        return "https://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/CAT_V2/assemblyHub/CHM13/CHM13.2bit"
    else:
        return f"https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome_build}/bigZips/{wildcards.genome_build}.2bit"


rule download_reference:
    output:
        ReferenceGenome.path(file_ending="/reference.2bit"),
    params:
        url = get_reference_url
    shell:
        "wget -O - {params.url} > {output}"


rule convert_reference_genome_to_fasta:
    input:
        ReferenceGenome.path(file_ending="/reference.2bit"),
    output:
        ReferenceGenome.path(),
    conda:
        "../envs/twobittofa.yml"
    params:
        remove_chr_prefix = lambda wildcards, input, output: " && sed -i 's/>chr/>/g' " + output[0] if config["genomes"][wildcards.genome_build]["remove_chr_prefix"] else ""
    shell:
        "twoBitToFa {input} {output} {params.remove_chr_prefix}"
    #wrapper:
        #    "v1.21.2/bio/ucsc/twoBitToFa"


rule get_chm13_reference:
    output:
        ReferenceGenome.path(genome_build="chm13"),
    shell:
        """
        cp local_data/CHM13v11Y.fa {output}
        """


ruleorder: get_chm13_reference > convert_reference_genome_to_fasta


rule get_dataset_reference:
    input:
        ref = ReferenceGenome.path(),
        index = ReferenceGenome.path(file_ending="/reference.fa.fai")
    output:
        tmp_genome = temp(BaseGenome.path(file_ending="/reference_tmp.fa")),
        #base_genome = BaseGenome.path()
    conda:
        "../envs/samtools.yml"
    params:
        regions=lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"].replace(",", " "),
    shell:
        "samtools faidx {input.ref} {params.regions} > {output.tmp_genome} "

rule process_dataset_reference:
    input:
        tmp_genome= BaseGenome.path(file_ending="/reference_tmp.fa"),
    output:
        base_genome = BaseGenome.path()
    conda:
        "../envs/python.yml"
    shell:
        "python scripts/format_fasta_headers.py {input.tmp_genome} {output.base_genome}"


# Simulates a set of variants that will be a source for a population
rule simulate_variant_source:
    input:
        base_genome = BaseGenome.path(),
        fai = BaseGenome.path(file_ending="/reference.fa.fai")
    output:
        variants = SimulatedVariantSource.path()
    params:
        tmp_output = lambda wildcards, input, output: output[0].replace(".vcf", ".tmp.vcf")
    #conda:
    #    "../envs/mason.yml"
    shell:
        "python scripts/simulate_variants.py {input.base_genome} "
        "{wildcards.snp_rate} "
        "{wildcards.small_indel_rate} "
        "{wildcards.sv_indel_rate} "
        "{output.variants}  "
        #"&& python3 scripts/remove_overlapping_indels.py {params.tmp_output} > {output}"


rule _outdated:
    shell:
        """
        mason_variator --seed 123 -ir {input.base_genome} -ov {params.tmp_output} \
        --snp-rate {wildcards.snp_rate} \
        --small-indel-rate {wildcards.small_indel_rate} \
        --sv-indel-rate {wildcards.sv_indel_rate} \
        --min-small-indel-size 1 \
        --max-small-indel-size 6 \
        --sv-inversion-rate 0 \
        --sv-translocation-rate 0 \
        --sv-duplication-rate 0 && \
        python3 scripts/remove_overlapping_indels.py {params.tmp_output} > {output}
        """

