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


rule pad_indels_for_paragraph:
    input:
        ref = BaseGenome.path(),
        vcf=FilteredPopulation.path(file_ending="/filtered_population_no_genotypes.vcf.gz"),
    output:
        vcf=FilteredPopulation.path(file_ending="/filtered_population_no_genotypes_padded.vcf"),
        vcf_gz=FilteredPopulation.path(file_ending="/filtered_population_no_genotypes_padded.vcf.gz")
    shell:
        """
        python scripts/fix_vcf_for_paragraph.py {input.ref} {input.vcf} {output.vcf} && bcftools sort {output.vcf} | bgzip -c > {output.vcf_gz}
        """



# Paragrap resorst to vcf with 1 individual, n_individuals doesn't matter for Paragraph, should give same output independelly of individuals in vcf
rule run_paragraph:
    input:
        sample_file=GenotypeResults.path(n_individuals=1, method="paragraph", file_ending="/samples.txt"),
        reference = BaseGenome.path(),
        variants = FilteredPopulation.path(n_individuals=1, file_ending="/filtered_population_no_genotypes_padded.vcf.gz"),
        index = FilteredPopulation.path(n_individuals=1, file_ending="/filtered_population_no_genotypes_padded.vcf.gz.tbi"),
        mapped_reads = Reads.path(file_ending="/reads.bam")
    output:
        results = GenotypeResults.path(n_individuals=1, method="paragraph", file_ending="/genotypes.vcf")
    params:
        output_folder = lambda wildcards, input, output: os.path.sep.join(output.results.split(os.path.sep)[:-1])
    threads:
        lambda wildcards: int(wildcards.n_threads)
    benchmark:
        GenotypeResults.path(method="paragraph", n_individuals=1, file_ending="/benchmark.csv")
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
        vcf=GenotypeResults.path(n_individuals=1,method="paragraph",file_ending="/genotypes.vcf"),
        benchmark=GenotypeResults.path(n_individuals=1,method="paragraph",file_ending="/benchmark.csv")
    output:
        vcf=GenotypeResults.path(method="paragraph",file_ending="/genotypes.vcf"),
        benchmark=GenotypeResults.path(method="paragraph",file_ending="/benchmark.csv")
    shell:
        "cp {input.vcf} {output.vcf} &&  cp {input.benchmark} {output.benchmark}"


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
