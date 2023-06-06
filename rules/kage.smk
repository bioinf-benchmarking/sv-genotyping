
rule remove_genotype_info:
    input: "{sample}.vcf"
    output: "{sample}_no_genotypes.vcf"
    shell: "cat {input} | cut -f 1-9 - > {output}"


rule kage_no_impuation_index:
    input:
        reference=BaseGenome.path(),
        population_vcf=PopulationWithoutIndividual.path(),
        population_vcf_no_genotypes=PopulationWithoutIndividual.path(file_ending="/population_no_genotypes.vcf"),
        reads=Reads.path(file_ending="/reads.fq")
    output:
        index = GenotypeResults.path(method="kage_no_imputation",file_ending="/index.npz")
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index} -V {input.population_vcf_no_genotypes} --modulo 200000033
        """


rule kage_index:
    input:
        reference = BaseGenome.path(),
        population_vcf = PopulationWithoutIndividual.path(),
        population_vcf_no_genotypes = PopulationWithoutIndividual.path(file_ending="/population_no_genotypes.vcf"),
        reads = Reads.path(file_ending="/reads.fq")
    output:
       index = GenotypeResults.path(method="kage", file_ending="/index.npz")
    shell:
        """
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index} -V {input.population_vcf_no_genotypes} --make-helper-model True --modulo 200000033
        """


rule run_kage_no_impuation:
    input:
        index = GenotypeResults.path(method="kage_no_imputation", file_ending="/index.npz"),
        reads = Reads.path()
    output:
        results = GenotypeResults.path(method="kage_no_imputation")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage}"


rule run_kage:
    input:
        index = GenotypeResults.path(method="kage", file_ending="/index.npz"),
        reads = Reads.path()
    output:
        results = GenotypeResults.path(method="kage")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --average-coverage {wildcards.coverage}"

