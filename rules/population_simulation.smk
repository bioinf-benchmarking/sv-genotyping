

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



# Removes an individual from a population
# Then also removes variants that only that individual had (using --min-ac 1)
rule remove_individual_from_population:
    input:
        population = RawPopulation.path(),
        sample_name = RawPopulation.path(file_ending="/unfiltered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        population = PopulationWithoutIndividual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        #"bcftools view --samples ^{wildcards.individual_id} -o {output.population} {input.population}"
        "bcftools view --samples-file ^{input.sample_name} {input.population} | bcftools view --min-ac 1 |"
        # add AF tag
        "bcftools +fill-tags -O z -o {output.population}  -- -t AF"


rule extract_individual_from_population:
    input:
        population = RawPopulation.path(),
        sample_name= RawPopulation.path(file_ending="/unfiltered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        individual = Individual.path()
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples-file {input.sample_name} {input.population} -o {output.individual}
        """
