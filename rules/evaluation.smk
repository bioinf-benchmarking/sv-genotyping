

def get_individual_for_genotype_accuracy(wildcards):
    print(Individual.from_flat_params(*wildcards).path(file_ending="/individual.vcf"))
    print(wildcards.method)
    print(Individual.path())
    print("NAMES")
    print([f.name for f in Individual.fields()])
    if "_multiallelic" in wildcards.method:
        return Individual.path(file_ending="/individual.multiallelic.vcf")
    else:
        return Individual.path()


def genotype_accuracy_input(wildcards):
    print(wildcards)
    if "pangenie" in wildcards.method and int(wildcards.n_individuals) > 125:
        return []
    else:
        out = [
            FilteredIndividual.path(file_ending="/individual.vcf"),
            StratifiedGenotypeResults.path()
        ]
        print(out)
        return out


rule genotype_accuracy:
    input:
        genotype_accuracy_input
        #truth = get_individual_for_genotype_accuracy,
        #truth = Individual.path(file_ending="/individual.vcf"),
        #genotypes = GenotypeResults.path()
    output:
        recall = GenotypeRecall.path(),
        one_minus_precision = GenotypeOneMinusPrecision.path(),
        #f1 = GenotypeF1Score.path(),
        report = GenotypeReport.path()
    script:
        "../scripts/genotype_accuracy.py"

rule f1_score:
    """ Using rtg tools"""
    input:
        rtg_report = VcfEvalReport.path()
    output:
        f1 = GenotypeF1Score.path(),
    run:
        with open(input.rtg_report) as f:
            line = f.readlines()[-1]
            f1_score = float(line.split()[-1])
            with open(output.f1, "w") as f:
                f.write(str(f1_score))

"""
rule genotype_accuracy_multiallelic:
    input:
        #truth = get_individual_for_genotype_accuracy,
        truth = Individual.path(file_ending="/individual.multiallelic.vcf"),
        genotypes = GenotypeResults.path(method="kage_multiallelic")
    output:
        f1 = GenotypeF1Score.path(method="kage_multiallelic"),
    script:
        "../scripts/genotype_accuracy.py"


ruleorder: genotype_accuracy_multiallelic > genotype_accuracy
"""

rule kage_debug:
    input:
        truth = Individual.path(),
        genotypes = GenotypeResults.path(),
        index=GenotypeResults.path(file_ending="/index.npz"),
        report = GenotypeReport.path(),
        node_counts = GenotypeResults.path(file_ending="/genotypes.vcf.node_counts.npy")
    output:
        debug = touch(GenotypeDebug.path())
    shell:
        "kage debug -i {input.index} -g {input.genotypes} -t {input.truth} -r {input.report} -n {input.node_counts}"




rule get_runtime:
    input:
        GenotypeResults.as_output(file_ending="/benchmark.csv")
        #f"data/{parameters.until('n_threads')}/benchmark.csv"
    output:
        GenotypeRuntime.path()
        #f"data/{parameters}/runtime.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 1 > {output}"


def get_stratification_file(wildcards):
    try:
        file = config["genome_stratification_files"][wildcards.genome_build][wildcards.stratification_type]
    except KeyError:
        print("Didn't find stratification file. Does the genome have stratification files in config.yaml?")
        raise
    return file


rule download_genome_stratification:
    output:
        GenomeStratification.path()
    params:
        file = get_stratification_file
    shell:
        "wget -O - {params.file} | gunzip -c > {output}"



rule subset_genotypes_on_stratification:
    input:
        vcf = GenotypeResults.path(),
        stratification=GenomeStratification.path()
    output:
        vcf = StratifiedGenotypeResults.path()
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools intersect -header -wa -u -a {input.vcf} -b {input.stratification} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} "


# "all" is simply the vcf copied
rule all_stratification:
    input:
        vcf = GenotypeResults.path(),
    output:
        vcf = StratifiedGenotypeResults.path(stratification_type="all")
    shell:
        "cat {input.vcf} | python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} "


ruleorder: all_stratification > subset_genotypes_on_stratification


rule subset_individual_on_stratification:
    input:
        vcf = Individual.path(),
        stratification=GenomeStratification.path()
    output:
        vcf = StratifiedIndividual.path()
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools intersect -header -wa -u -a {input.vcf} -b {input.stratification} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} "


# "all" is simply the vcf copied
rule all_stratification_individual:
    input:
        vcf = Individual.path(),
    output:
        vcf = StratifiedIndividual.path(stratification_type="all")
    shell:
        "cat {input.vcf} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} "


ruleorder: all_stratification_individual > subset_individual_on_stratification



rule filter_individual:
    """
    Can filter an individual so that it only keeps variants that are in the population
    """
    input:
        individual = StratifiedIndividual.path(),
        population = PopulationWithoutIndividual.path(),
    output:
        individual = FilteredIndividual.path(individual_filter="only_variants_in_population")
    params:
        dir = lambda wildcards, input, output: "/".join(output.individual.split("/")[:-1])
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bgzip -c {input.individual} > {input.individual}.gz &&
        tabix -p vcf {input.individual}.gz &&
        bcftools isec -p {params.dir} {input.individual}.gz {input.population} -O z && 
        zcat {params.dir}/0002.vcf.gz > {output.individual}
        """


rule filter_individual_no_filter:
    input:
        individual = StratifiedIndividual.path(),
    output:
        individual = FilteredIndividual.path(individual_filter="none")
    shell:
        "cp {input.individual} {output.individual}"



rule rtg_format:
    input:
        BaseGenome.path()
    output:
        BaseGenome.path(file_ending="/sdf/done")
    conda:
        "../envs/rtg_tools.yml"
    params:
        out_path = lambda wildcards, input, output: "/".join(output[0].split("/")[:-1])
    shell:
        """
        rtg format -o {params.out_path} {input}
        """


rule rtg_tools_vcf_eval:
    input:
        genotypes=FilteredIndividual.path(file_ending="/individual.vcf"),
        truth=StratifiedGenotypeResults.path(),
        sdf=BaseGenome.path(file_ending="/sdf/done")
    output:
        report = VcfEvalReport.path()
    params:
        out_dir = lambda wildcards, input, output: "/".join(output.report.split("/")[:-1]) + "/rtg_tools_tmp",
        sdf = lambda wildcards, input, output: input.sdf.replace("/done", "")
    conda:
        "../envs/rtg_tools.yml"
    shell:
        """
        bgzip -c {input.genotypes} > {input.genotypes}.gz &&
        tabix -p vcf {input.genotypes}.gz &&
        bgzip -c {input.truth} > {input.truth}.gz &&
        tabix -p vcf {input.truth}.gz &&
        rm -rf {params.out_dir} &&
        
        rtg vcfeval -b {input.truth}.gz -c {input.genotypes}.gz -t {params.sdf} -o {params.out_dir} --no-roc &&
        mv {params.out_dir}/summary.txt {output.report}
        """
