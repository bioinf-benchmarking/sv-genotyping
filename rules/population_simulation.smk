

rule simulate_population:
    input:
        variant_source = SimulatedVariantSource.path()
    output:
        population = SimulatedPopulation.path()
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
        population = Population.path(),
        sample_name = Population.path(file_ending="/population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        population = PopulationWithoutIndividual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        #"bcftools view --samples ^{wildcards.individual_id} -o {output.population} {input.population}"
        "bcftools view --samples-file ^{input.sample_name} -o {output.population} {input.population}"


rule extract_individual_from_population:
    input:
        population = Population.path(),
        sample_name= Population.path(file_ending="/population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        individual = Individual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples-file {input.sample_name} -o {output.individual} {input.population}
        """
