
def get_chromosomes(wildcards):
    out = ",".join([c.split(":")[0] for c in config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"].split(",")])
    return out


rule run_glimpse:
    input:
        vcf = GenotypeResults.path(method="kage_no_imputation"),
        ref_vcf = FilteredPopulation.path(),
    output:
        vcf = GenotypeResults.path(method="kage_with_glimpse"),
    benchmark:
        GenotypeResults.path(method="kage_with_glimpse", file_ending="/glimpse.csv")
    threads: lambda wildcards: int(wildcards.n_threads)
    params:
        regions = get_chromosomes,
    shell:
        """
        kage glimpse -p {input.ref_vcf} -g {input.vcf} -o {output.vcf} -t {threads} -c {params.regions}
        """