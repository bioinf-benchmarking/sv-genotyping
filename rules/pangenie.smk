# Modified version of Pangenie run-from-callset snakemake pipeline (https://bitbucket.org/jana_ebler/pangenie/raw/d340042099a3af7683729327235a3cd3c8bb8a71/pipelines/run-from-callset/Snakefile)

# assign IDs to all alleles
rule add_ids:
    input:
        vcf="{file}.vcf"
        #vcf=PopulationWithoutIndividual.path()
    output:
        vcf = "{file}.with_ids.vcf"
        #vcf=PopulationWithoutIndividual.path(file_ending="/population.with_ids.vcf"),
    shell:
        'cat {input} | python3 scripts/pangenie_add_ids.py > {output}'


# merge variants into a pangenome graph
rule merge_haplotypes:
    input:
        vcf=PopulationWithoutIndividual.path(file_ending="/population.with_ids.vcf"),
        reference = BaseGenome.path(),
    output:
        population_vcf = PopulationWithoutIndividual.path(file_ending="/population.multiallelic.vcf"),
    shell:
        """
        python3 scripts/pangenie_merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 > {output}
        """


rule merge_individual_haplotypes:
    input:
        vcf=Individual.path(file_ending="/individual.with_ids.vcf"),
        reference = BaseGenome.path(),
    output:
        vcf = Individual.path(file_ending="/individual.multiallelic.vcf"),
    shell:
        """
        python3 scripts/pangenie_merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 > {output}
        """


rule run_pangenie:
    input:
        reference = BaseGenome.path(),
        population_vcf = PopulationWithoutIndividual.path(file_ending="/population.multiallelic.vcf"),
        reads = Reads.path(file_ending="/reads.fq.gz")
    output:
        results = GenotypeResults.path(method="pangenie", file_ending="/genotypes_multiallelic.vcf")
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


# rule that outputs only the multiallelic vcf from pangenie as a "pangenie_multiallelic" method
rule run_pangenie_multiallelic:
    input:
        results = GenotypeResults.path(method="pangenie",file_ending="/genotypes_multiallelic.vcf")
    output:
        results = GenotypeResults.path(method="pangenie_multiallelic",file_ending="/genotypes.vcf")
    shell:
        "cp {input} {output}"


rule convert_pangenie_genotypes_back_to_biallelic:
    input:
        callset_vcf=PopulationWithoutIndividual.path(file_ending="/population.with_ids.vcf"),
        vcf = GenotypeResults.path(method="pangenie",file_ending="/genotypes_multiallelic.vcf")
    output:
        vcf = GenotypeResults.path(method="pangenie",file_ending="/genotypes.vcf")
    shell:
        """
        cat {input.vcf} | python scripts/pangenie_convert_to_biallelic.py {input.callset_vcf} > {output.vcf}
        """



