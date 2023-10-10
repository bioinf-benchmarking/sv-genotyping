import os

rule make_paragraph_samples_file:
    input:
        sorted_bam=Reads.path(file_ending="/reads.sorted.bam")
    output:
        #sample_file="data/{dataset}/paragraph_samples_{sample_id,[a-zA-Z0-9]+}_{read_info}_{coverage,\d+}x.txt"
        sample_file=GenotypeResults.path(method="paragraph", file_ending="/samples.txt")
    shell:
        """
        echo 'id\tpath\tdepth\tread length\n{wildcards.individual_id}\t{input.sorted_bam}\t{wildcards.coverage}\t{wildcards.read_length}' > {output.sample_file}
        """


# Paragrap resorst to vcf with 1 individual, n_individuals doesn't matter for Paragraph, should give same output independelly of individuals in vcf
rule run_paragraph:
    input:
        sample_file=GenotypeResults.path(n_individuals=1, method="paragraph", file_ending="/samples.txt"),
        reference = BaseGenome.path(),
        variants = FilteredPopulation.path(n_individuals=1, file_ending="/filtered_population_no_genotypes.vcf.gz"),
        index = FilteredPopulation.path(n_individuals=1, file_ending="/filtered_population_no_genotypes.vcf.gz.tbi"),
        mapped_reads = Reads.path(file_ending="/reads.bam")
    output:
        results = GenotypeResults.path(n_individuals=1, method="paragraph", file_ending="/genotypes.vcf")
    params:
        output_folder = lambda wildcards, input, output: os.path.sep.join(output.results.split(os.path.sep)[:-1])
    threads:
        lambda wildcards: int(wildcards.n_threads)
    conda:
        "../envs/paragraph.yml"
    shell:
        """
        python paragraph/bin/multigrmpy.py -i {input.variants} \
        -m {input.sample_file} \
        -r {input.reference} \
        -o {params.output_folder} \
        -t {config[n_threads]} &&
        zcat {params.output_folder}/genotypes.vcf.gz > {output.results}
        """


rule run_paragraph_wrapper:
    input:
        GenotypeResults.path(n_individuals=1,method="paragraph",file_ending="/genotypes.vcf")
    output:
        GenotypeResults.path(method="paragraph",file_ending="/genotypes.vcf")
    shell:
        "cp {input} {output} "


"""
rule convert_paragraph_genotypes_back_to_biallelic:
    input:
        callset_vcf=FilteredPopulation.path(file_ending="/filtered_population.with_ids.vcf"),
        vcf = GenotypeResults.path(method="paragraph",file_ending="/genotypes_multiallelic.vcf")
    output:
        vcf = GenotypeResults.path(method="paragraph",file_ending="/genotypes.vcf")
    shell:
        "cat {input.vcf} | python scripts/pangenie_convert_to_biallelic.py {input.callset_vcf} > {output.vcf}"
"""
