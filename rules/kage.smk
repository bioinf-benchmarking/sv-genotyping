
rule remove_genotype_info:
    input: "{sample}.vcf"
    output: "{sample}_no_genotypes.vcf"
    shell: "cat {input} | cut -f 1-9 - > {output}"


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
        kage index -r {input.reference} -v {input.population_vcf} -o {output.index} -V {input.population_vcf_no_genotypes}
        """


rule run_kage:
    input:
        index = GenotypeResults.path(method="kage", file_ending="/index.npz"),
        reads = Reads.path()
    output:
        results = GenotypeResults.path(method="kage")
    shell:
        "kage genotype -i {input.index} -r {input.reads} -o {output.results} -t {wildcards.n_threads} --ignore-helper-model True --ignore-helper-variants True"