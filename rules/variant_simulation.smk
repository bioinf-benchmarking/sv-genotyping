
rule test:
    output:
        touch("test.txt")



rule get_dataset_reference:
    input:
        ref = ReferenceGenome.path(),
        index = ReferenceGenome.path(file_ending="/reference.fa.fai")
    output:
        base_genome = BaseGenome.path()
    conda:
        "../envs/samtools.yml"
    params:
        regions=lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"].replace(",", " ")
    shell:
        "samtools faidx {input.ref} {params.regions} > {output.base_genome}"


# Simulates a set of variants that will be a source for a population
rule simulate_variant_source:
    input:
        base_genome = BaseGenome.path()
    output:
        variants = VariantSource.path()
    params:
        tmp_output = lambda wildcards, input, output: output[0].replace(".vcf", ".tmp.vcf")
    conda:
        "../envs/mason.yml"
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

