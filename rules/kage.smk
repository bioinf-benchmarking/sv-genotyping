
rule remove_genotype_info:
    input: "{sample}.vcf"
    output: "{sample}_no_genotypes.vcf"
    shell: "cat {input} | cut -f 1-9 -d$'\t' - > {output}"


rule remove_genotype_info_gz:
    input: "{sample}.vcf.gz"
    output: "{sample}_no_genotypes.vcf.gz"
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view -G {input} -O z -o {output}"



rule kage_index:
    input:
        reference = BaseGenome.path(),
        population_vcf = FilteredPopulation.path(),
    output:
       index = GenotypeResults.path(method="kage", file_ending="/index.npz")
    threads:
        lambda wildcards: int(wildcards.n_threads)
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index}  \
        --make-helper-model True --modulo 200000033 --variant-window 5 -k 31
        """


rule run_kage_no_impuation:
    input:
        index = GenotypeResults.path(method="kage", file_ending="/index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage_no_imputation")
    benchmark:
        GenotypeResults.path(method="kage_no_imputation", file_ending="/benchmark.csv")
    shell:
        "kage genotype -i {input.index} -r {input.reads} "
        " -o {output.results} -t {wildcards.n_threads} "
        "--average-coverage {wildcards.coverage} -k 31 "
        "--ignore-helper-model True "


rule run_kage:
    input:
        index = GenotypeResults.path(method="kage", file_ending="/index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage"),
        node_counts = GenotypeResults.path(method="kage", file_ending="/genotypes.vcf.node_counts.npy")
    benchmark:
        GenotypeResults.path(method="kage", file_ending="/benchmark.csv")
    threads:
        lambda wildcards: int(wildcards.n_threads)
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage} -k 31"

