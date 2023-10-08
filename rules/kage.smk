
rule remove_genotype_info:
    input: "{sample}.vcf"
    output: "{sample}_no_genotypes.vcf"
    shell: "cat {input} | cut -f 1-9 - > {output}"


rule kage_no_impuation_index:
    input:
        reference=BaseGenome.path(),
        population_vcf=PopulationWithoutIndividual.path(),
        population_vcf_no_genotypes=PopulationWithoutIndividual.path(file_ending="/population_no_genotypes.vcf"),
    output:
        index = GenotypeResults.path(method="kage_no_imputation",file_ending="/index.npz")
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index} -V {input.population_vcf_no_genotypes} --modulo 200000033 --variant-window 4
        """


rule kage_index:
    input:
        reference = BaseGenome.path(),
        population_vcf = PopulationWithoutIndividual.path(),
    output:
       index = GenotypeResults.path(method="kage", file_ending="/index.npz")
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index}  --make-helper-model True --modulo 200000033 --variant-window 7 -k 31
        """


# Runs kage index on a multiallelic vcf (for testing that that also works)
rule kage_index_multiallelic:
    input:
        reference = BaseGenome.path(),
        population_vcf = PopulationWithoutIndividual.path(file_ending="/population.multiallelic.vcf"),
    output:
       index = GenotypeResults.path(method="kage_multiallelic", file_ending="/index.npz")
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index}  --make-helper-model True --modulo 200000033 --variant-window 6 -k 31
        """


rule run_kage_no_impuation:
    input:
        index = GenotypeResults.path(method="kage_no_imputation", file_ending="/index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage_no_imputation")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage}"


rule run_kage:
    input:
        index = GenotypeResults.path(method="kage", file_ending="/index.npz"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage"),
        node_counts = GenotypeResults.path(method="kage", file_ending="/genotypes.vcf.node_counts.npy")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage} -k 31"



rule run_kage_multiallelic:
    input:
        index=GenotypeResults.path(method="kage_multiallelic", file_ending="/index.npz"),
        reads=Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="kage_multiallelic"),
        node_counts=GenotypeResults.path(method="kage_multiallelic", file_ending="/genotypes.vcf.node_counts.npy")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage} -k 31"

