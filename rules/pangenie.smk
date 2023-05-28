

rule run_pangenie:
    input:
        reference = BaseGenome.path(),
        population_vcf = PopulationWithoutIndividual.path(),
        reads = Reads.path(file_ending="/reads.fq")
    output:
        results = GenotypeResults.path(method="pangenie")
    params:
        jellyfish_memory = lambda wildcards: 3000000000 if wildcards.size == "big" else 300000000
    shell:
        """
        PanGenie -i {input.reads} \
        -r {input.reference} \
        -v {input.population_vcf} \
        -e {params.jellyfish_memory} \
        -o {output.results} \
        -t {wildcards.n_threads} \
        -j {wildcards.n_threads} \
        && mv {output.results}_genotyping.vcf {output.results}
        """