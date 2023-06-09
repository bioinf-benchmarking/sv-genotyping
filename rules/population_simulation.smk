

rule simulate_population:
    input:
        variant_source = VariantSource.path()
    output:
        population = Population.path()
    shell:
        """
        python3 scripts/population_simulation.py {input.variant_source} \
        {output.population} \
        {wildcards.n_individuals} \
        {wildcards.allele_frequency} \
        {wildcards.correlation} 
        """


rule remove_individual_from_population:
    input:
        population = Population.path()
    output:
        population = PopulationWithoutIndividual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples ^{wildcards.individual_id} -o {output.population} {input.population}
        """


rule extract_individual_from_population:
    input:
        population=Population.path()
    output:
        individual = Individual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples {wildcards.individual_id} -o {output.individual} {input.population}
        """
