

rule simulate_population:
    input:
        variant_source = SimulatedVariantSource.path()
    output:
        tmp =  SimulatedPopulation.path(file_ending="/unfiltered_population.tmp.vcf.gz"),
        population = SimulatedPopulation.path()
    shell:
        """
        python3 scripts/population_simulation.py {input.variant_source} \
        {output.tmp} \
        3000 \
        {wildcards.allele_frequency} \
        {wildcards.correlation} \
        && bgzip -c {output.tmp} > {output.population}
        """



rule remove_individual_from_population:
    input:
        population = FilteredPopulation.path(),
        sample_name = FilteredPopulation.path(file_ending="/filtered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        population = PopulationWithoutIndividual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        #"bcftools view --samples ^{wildcards.individual_id} -o {output.population} {input.population}"
        "bcftools view --samples-file ^{input.sample_name} -o {output.population} {input.population}"


rule extract_individual_from_population:
    input:
        population = FilteredPopulation.path(),
        sample_name= FilteredPopulation.path(file_ending="/filtered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        individual = Individual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples-file {input.sample_name} -o {output.individual} {input.population}
        """
