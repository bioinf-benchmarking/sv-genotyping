# Modified version of Pangenie run-from-callset snakemake pipeline (https://bitbucket.org/jana_ebler/pangenie/raw/d340042099a3af7683729327235a3cd3c8bb8a71/pipelines/run-from-callset/Snakefile)

# assign IDs to all alleles
rule add_ids:
    input:
        vcf="{file}.vcf.gz"
        #vcf=PopulationWithoutIndividual.path()
    output:
        vcf = "{file}.with_ids.vcf"
        #vcf=PopulationWithoutIndividual.path(file_ending="/population.with_ids.vcf"),
    shell:
        'zcat {input} | python3 scripts/pangenie_add_ids.py > {output}'


# merge variants into a pangenome graph
rule merge_haplotypes:
    input:
        vcf=FilteredPopulation.path(file_ending="/filtered_population.with_ids.vcf"),
        reference = BaseGenome.path(),
    output:
        population_vcf = FilteredPopulation.path(file_ending="/filtered_population.multiallelic.vcf"),
        population_vcfgz = FilteredPopulation.path(file_ending="/filtered_population.multiallelic.vcf.gz"),
    shell:
        """
        python3 scripts/pangenie_merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 > {output.population_vcf} &&
        bgzip -c {output.population_vcf} > {output.population_vcfgz} 
        """


rule merge_individual_haplotypes:
    input:
        vcf=Individual.path(file_ending="/individual.with_ids.vcf"),
        reference = BaseGenome.path(),
    output:
        vcf = Individual.path(file_ending="/individual.multiallelic.vcf"),
        vcfgz= Individual.path(file_ending="/individual.multiallelic.vcf.gz"),
    shell:
        """
        python3 scripts/pangenie_merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 > {output} &&
        bgzip -c {output.vcf} > {output.vcfgz} 
        """


rule pangenie_index:
    input:
        reference = BaseGenome.path(),
        population_vcf = FilteredPopulation.path(file_ending="/filtered_population.multiallelic.vcf"),
    output:
        FilteredPopulation.path(file_ending=["/pangenie_index_path_segments.fasta", "/pangenie_index_UniqueKmersMap.cereal"])
    params:
        out_prefix = lambda wildcards, input, output: output[0].replace("index_path_segments.fasta", "index")
    threads:
        8
    shell:
        """
        PanGenie-index -v {input.population_vcf} -r {input.reference} -t 8 -o {params.out_prefix}
        """


rule run_pangenie:
    input:
        #reference = BaseGenome.path(),
        #population_vcf = FilteredPopulation.path(file_ending="/filtered_population.multiallelic.vcf"),
        index = FilteredPopulation.path(file_ending="/pangenie_index_path_segments.fasta"),
        reads = Reads.path(file_ending="/reads.fq")
    output:
        results = GenotypeResults.path(method="pangenie", file_ending="/genotypes_multiallelic.vcf")
    params:
        index_prefix = lambda wildcards, input, output: input.index.replace("index_path_segments.fasta", "index"),
        jellyfish_memory = lambda wildcards: 3000000000 if wildcards.size == "big" else 300000000
    threads:
        lambda wildcards: int(wildcards.n_threads)
    benchmark:
        GenotypeResults.path(method="pangenie", file_ending="/benchmark.csv")
    shell:
        """
         PanGenie -i {input.reads} \
         -f {params.index_prefix} \
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
        callset_vcf=FilteredPopulation.path(file_ending="/filtered_population.with_ids.vcf"),
        vcf = GenotypeResults.path(method="pangenie",file_ending="/genotypes_multiallelic.vcf")
    output:
        vcf = GenotypeResults.path(method="pangenie",file_ending="/genotypes.vcf")
    shell:
        """
        cat {input.vcf} | python scripts/pangenie_convert_to_biallelic.py {input.callset_vcf} > {output.vcf}
        """



