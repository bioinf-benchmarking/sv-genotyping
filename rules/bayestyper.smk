
rule run_kmc_bayestyper:
    input:
        reads = Reads.path()
    output:
        pre=Reads.path(file_ending="/kmers_bayestyper.kmc_pre"),
        suf=Reads.path(file_ending="/kmers_bayestyper.kmc_suf")
    threads:
        config["n_threads"]
    resources:
        mem_gb=20
    benchmark:
        Reads.path(file_ending="/kmc.csv")
        #"data/{dataset}/benchmarks/bayestyper_kmc_{reads}.tsv"
    params:
        out_prefix = lambda wildcards, input, output: output[0].replace(".kmc_pre", "")
    conda:
        "../envs/kmc.yml"
    shell:
        "mkdir -p {params.out_prefix}/kmc_tmp_bayestyper && "
        "kmc -t{config[n_threads]} -k55 -ci1 -fq {input} {params.out_prefix} {params.out_prefix}/kmc_tmp_bayestyper"


rule make_multiallelic_variants_for_bayestyper:
    input:
        vcf=FilteredPopulationWithCollapsedAlleles.path(file_ending="/filtered_population_with_collapsed_alleles_no_genotypes.vcf.gz"),
        ref = BaseGenome.path()
    output:
        FilteredPopulationWithCollapsedAlleles.path(file_ending="/filtered_population_with_collapsed_alleles_no_genotypes_multiallelic.vcf"),
    conda: "../envs/bcftools118.yml"
    shell:
        # bcftools norm twice, likely something buggy in bcftools, doesn't sort properly the first time and misses some
        "bcftools norm -m +any -f {input.ref} {input.vcf} | bcftools sort | bcftools norm -m +any -f {input.ref} > {output}"


rule make_samples_tsv_for_bayestyper:
    input:
        pre = Reads.path(file_ending="/kmers_bayestyper.kmc_pre"),
    output:
        samples=Reads.path(file_ending="/samples.tsv"),
        #"data/{dataset}/samples_{individual,\w+}_{reads_type}_{coverage,\d+}x.tsv"
    params:
        kmers_prefix = lambda wildcards, input, output: input[0].replace(".kmc_pre", "")
    shell:
        "echo '{wildcards.individual_id}\tM\t{params.kmers_prefix}' > {output}"


rule make_bloomfilter_for_bayestyper:
    input:
        pre=Reads.path(file_ending="/kmers_bayestyper.kmc_pre"),
        suf=Reads.path(file_ending="/kmers_bayestyper.kmc_suf")
    output:
        bloomdata=Reads.path(file_ending="/kmers_bayestyper.bloomData"),
        bloommeta=Reads.path(file_ending="/kmers_bayestyper.bloomMeta")
        #"data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        #"data/{dataset}/{reads}.kmers_bayestyper.bloomMeta"
    threads:
        config["n_threads"]
    benchmark:
        Reads.path(file_ending="/bloomfilter.csv")
        #"data/{dataset}/benchmarks/bayestyper_bloomfilter_{reads}.tsv"
    params:
        kmers_prefix = lambda wildcards, input, output: input[0].replace(".kmc_pre", "")
    conda: "../envs/bayestyper.yml"
    shell:
        "bayesTyperTools makeBloom -k {params.kmers_prefix} -p {config[n_threads]}"


rule make_decoy_fasta:
    output: BaseGenome.path(file_ending="/decoy.fasta")
    shell: "echo -n '' > {output}"

rule run_bayestyper:
    input:
        ref=BaseGenome.path(),  # "data/{dataset}/ref.fa",
        decoy=BaseGenome.path(file_ending="/decoy.fasta"),  #"data/{dataset}/decoy.fasta",
        bloomdata=Reads.path(file_ending="/kmers_bayestyper.bloomData"),
        bloommeta=Reads.path(file_ending="/kmers_bayestyper.bloomMeta"),
        #bloomdata="data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        #bloommeta="data/{dataset}/{reads}.kmers_bayestyper.bloomMeta",
        variants=FilteredPopulationWithCollapsedAlleles.path(file_ending="/filtered_population_with_collapsed_alleles_no_genotypes_multiallelic.vcf"),
        #variants="data/{dataset}/variants_no_genotypes_multiallelic.vcf",
        samples=Reads.path(file_ending="/samples.tsv"),
        #samples="data/{dataset}/samples_{reads}.tsv",
    output:
        #units=dynamic("data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_{unit_id}/variant_clusters.bin")
        #units="data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_1/variant_clusters.bin"
        #genotypes="data/{dataset,\w+}/bayestyper_{reads}.vcf.gz"
        genotypes=GenotypeResults.path(method="bayestyper", file_ending="/genotypes_multiallelic.vcf.gz")
    benchmark:
        GenotypeResults.path(method="bayestyper", file_ending="/running_bayestyper.csv")
    params:
        out_prefix="bayestyper/bayestyper",
        tmp_dir=lambda wildcards, input, output: "/".join(output[0].split("/")[:-1]) + "/tmp/"
    threads:
        lambda wildcards: int(wildcards.n_threads)
    resources:
        mem_gb=50
    conda: "../envs/bayestyper.yml"
    shell:
        "mkdir -p {params.tmp_dir} && "
        "rm -rf {params.tmp_dir}/* && "
        "bayesTyper cluster -v {input.variants} -s {input.samples} -g {input.ref} -d {input.decoy} -p {wildcards.n_threads} -o {params.tmp_dir}/bayestyper && "
        "for dir in {params.tmp_dir}/bayestyper_unit_*; do\n "
        "    bayesTyper genotype -v $dir/variant_clusters.bin -c {params.tmp_dir}/bayestyper_cluster_data -s {input.samples} -g {input.ref} -d {input.decoy} -o $dir/bayestyper -z -p {wildcards.n_threads}; "
        "done && "
        # hack since bcftools seem to not bgzip, but only gzip:
        "ls {params.tmp_dir}/bayestyper_unit_*/*.vcf.gz | xargs -P 16 -n 1 gunzip && "
        "ls {params.tmp_dir}/bayestyper_unit_*/*.vcf | xargs -P 16 -n 1 bgzip && "
        "ls {params.tmp_dir}/bayestyper_unit_*/*.vcf.gz | xargs -P 16 -n 1 tabix -f -p vcf && "

        "bcftools concat -O z -a -o {output.genotypes} {params.tmp_dir}/bayestyper_unit_*/bayestyper.vcf.gz"


rule convert_bayestyper_multiallelic_to_biallelic:
    input:
        genotypes=GenotypeResults.path(method="bayestyper", file_ending="/genotypes_multiallelic.vcf.gz"),
        ref=BaseGenome.path()
    output:
        genotypes=GenotypeResults.path(method="bayestyper", file_ending="/genotypes.vcf")
    conda: "../envs/bcftools.yml"
    shell:
        "bcftools norm -m -any -f {input.ref} {input.genotypes} > {output.genotypes}"


"""
rule run_bayestyper_wrapper:
    input:
        vcf=GenotypeResults.path(n_individuals=1,method="bayestyper",file_ending="/genotypes.vcf"),
    output:
        vcf=GenotypeResults.path(method="bayestyper",file_ending="/genotypes.vcf"),
    shell:
        "cp {input.vcf} {output.vcf}"
"""

rule add_bayestyper_benchmark_times:
    input:
        bloomfilter=Reads.path(file_ending="/bloomfilter.csv"),
        kmc=Reads.path(file_ending="/kmc.csv"),
        bayestyper=GenotypeResults.path(method="bayestyper", file_ending="/running_bayestyper.csv")
    output:
        GenotypeRuntime.path(method="bayestyper")
    run:
        with open(output[0], "w") as f:
            total_time = 0
            for file_name in input:
                with open(file_name) as file:
                    time = float(file.readlines()[-1].split()[0])
                    print(file_name, time)
                    total_time += time
            f.write(str(total_time)+"\n")
