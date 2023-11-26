
rule remove_genotype_info:
    input: "{sample}.vcf"
    output: "{sample}_no_genotypes.vcf"
    shell: "cat {input} | cut -f 1-9 -d$'\t' - > {output}"


rule remove_genotype_info_gz:
    input:
        vcf="{sample}.vcf.gz",
        index="{sample}.vcf.gz.tbi"
    output:
        "{sample}_no_genotypes.vcf.gz"
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view -G {input.vcf} -O z -o {output}"



rule kage_index:
    input:
        reference = BaseGenome.path(),
        population_vcf = FilteredPopulationWithCollapsedAlleles.path(),
        population_vcf_no_genotypes = FilteredPopulationWithCollapsedAlleles.path(file_ending="/filtered_population_with_collapsed_alleles_no_genotypes.vcf.gz"),
    output:
       index = FilteredPopulationWithCollapsedAlleles.path(file_ending="/kage_index.npz")
    threads:
        lambda wildcards: 8
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -V {input.population_vcf_no_genotypes} -o {output.index}  \
        --modulo 200000033 --variant-window 7 -k 31 --min-af-deletions-filter 0.1
        """


rule run_kage_no_imputation:
    input:
        index = FilteredPopulationWithCollapsedAlleles.path(file_ending="/kage_index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage_no_imputation"),
        node_counts= GenotypeResults.path(method="kage_no_imputation",file_ending="/genotypes.vcf.node_counts.npy")
    benchmark:
        GenotypeResults.path(method="kage_no_imputation", file_ending="/benchmark.csv")
    shell:
        "kage genotype -i {input.index} -r {input.reads} "
        " -o {output.results} -t {wildcards.n_threads} "
        "--average-coverage {wildcards.coverage} -k 31 "
        "--ignore-helper-model True --write-debug-data True"


rule run_kage:
    input:
        index = FilteredPopulationWithCollapsedAlleles.path(file_ending="/kage_index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage"),
        node_counts = GenotypeResults.path(method="kage", file_ending="/genotypes.vcf.node_counts.npy")
    benchmark:
        GenotypeResults.path(method="kage", file_ending="/benchmark.csv")
    threads:
        lambda wildcards: int(wildcards.n_threads)
    shell:
        "kage genotype -i {input.index} "
        "-r {input.reads} "
        "-o {output.results} "
        "-t {wildcards.n_threads} "
        "--average-coverage {wildcards.coverage} "
        "-k 31 --write-debug-data True"
        #"--ignore-homo-ref True"  # don't write homo ref variants for faster writing

