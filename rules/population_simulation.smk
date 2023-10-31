

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
        population = PopulationWithoutIndividual.path(individual_source="from_pangenome")
    conda:
        "../envs/bcftools.yml"
    shell:
        #"bcftools view --samples ^{wildcards.individual_id} -o {output.population} {input.population}"
        "bcftools view --samples-file ^{input.sample_name} {input.population} | bcftools view --min-ac 1 |"
        # add AF tag
        "bcftools +fill-tags -O z -o {output.population}  -- -t AF"


rule remove_remote_individual_from_population:
    """
    Removes an invididual from population when individual is "remote", i.e. not specified by a random id in the population
    In these cases, it is assuemd that the individual is not in the population at all.
    """
    input:
        population=RawPopulation.path(),
    output:
        population = PopulationWithoutIndividual.path(individual_source="remote")
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        cp {input.population} {output.population} &&
        bcftools +fill-tags -O z -o {output.population} {input.population} -- -t AF
        """



rule extract_individual_from_population:
    input:
        population = RawPopulation.path(),
        sample_name= RawPopulation.path(file_ending="/unfiltered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        individual = Individual.path(individual_source="from_pangenome")
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --samples-file {input.sample_name} {input.population} -o {output.individual}
        """


rule download_individual:
    output:
        individual = Individual.path(individual_source="remote")
    conda:
        "../envs/bcftools.yml"
    params:
        chromosomes= lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"],
        url = lambda wildcards: config["real_individuals"][wildcards.individual_id]
    shell:
        """
        wget -O - {params.url} | zcat > {output.individual}.tmp.vcf &&
        sed 's/chr//g' {output.individual}.tmp.vcf | bgzip -c > {output.individual}.tmp.vcf.gz &&
        tabix -f -p vcf {output.individual}.tmp.vcf.gz &&
        bcftools view --regions {params.chromosomes} {output.individual}.tmp.vcf.gz > {output.individual}
        """
