
rule download_polaris_svs:
    output:
        RealVariantSource.path(database_name="polaris", variant_type="svs", file_ending="/raw.vcf.gz")
    shell:
        "zcat local_data/polaris.vcf.gz | sed 's/chr//g' | bgzip -c > {output}"


#only keeping SVs with known sequences
rule filter_polaris_svs:
    input:
        variants = RealVariantSource.path(database_name="polaris", variant_type="svs", file_ending="/raw.vcf.gz"),
        reference = ReferenceGenome.path()
    output:
        RealVariantSource.path(database_name="polaris", variant_type="svs", file_ending="/raw_without_unknown_sequences.vcf.gz")
    #conda:
    #    "../envs/bcftools.yml"
    shell:
        #"python scripts/filter_svs.py {input} | bgzip -c > {output}"
        "kage preprocess_sv_vcf -v {input.variants} -f {input.reference} | bgzip -c > {output}"

