

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
            BestGenotypes.path()
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
        f1 = GenotypeF1Score.path(),
        weighted_genotype_concordance = WeightedGenotypeConcordance.path(),
        report = GenotypeReport.path()
    script:
        "../scripts/genotype_accuracy.py"

"""
rule f1_score:
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
        index = FilteredPopulationWithCollapsedAlleles.path(file_ending="/kage_index.npz"),
        report = GenotypeReport.path(),
        node_counts = GenotypeResults.path(file_ending="/genotypes.vcf.node_counts.npy")
    output:
        debug = touch(GenotypeDebug.path())
    shell:
        "kage debug -i {input.index} -g {input.genotypes} -t {input.truth} -r {input.report} -n {input.node_counts} --probs {input.genotypes}.probs.npy --count-probs {input.genotypes}.count_probs.npy --numeric-genotypes {input.genotypes}.genotypes.npy"




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
        file = get_stratification_file,
        remove_chr_prefix= lambda wildcards,input,output: " && sed -i 's/chr//g' " + output[0] if config["genomes"][wildcards.genome_build]["remove_chr_prefix"] else "",
        gunzip = lambda wildcards: " | gunzip -c" if get_stratification_file(wildcards).endswith(".gz") else ""
    shell:
        "wget -O - {params.file} {params.gunzip} > {output} {params.remove_chr_prefix}"


rule download_syndip_genome_stratification:
    input:
        "data/syndip.tar"
    output:
        GenomeStratification.path(stratification_type="syndip-confident-regions")
    params:
        remove_chr_prefix= lambda wildcards,input,output: " && sed -i 's/chr//g' " + output[0] if config["genomes"][wildcards.genome_build]["remove_chr_prefix"] else "",
    shell:
        "tar -xvf syndip.tar && "
        "zcat CHM-eval.kit/full.38.bed.gz > {output} {params.remove_chr_prefix}"


ruleorder: download_syndip_genome_stratification > download_genome_stratification


rule subset_genome_stratification_on_dataset:
    input:
        stratification=GenomeStratification.path(),
        reference=BaseGenome.path(file_ending="/reference.fa.fai")
    output:
        stratification=GenomeStratificationOnDataset.path()
    run:
        with open(input.reference) as ref:
            chromosomes = set([line.split()[0] for line in ref.readlines()])

        with open(input.stratification) as strat:
            with open(output.stratification, "w") as out:
                for line in strat.readlines():
                    if line.split()[0] in chromosomes:
                        out.write(line)


# the "all" stratification, should just create a file of the whole dataset region
rule subset_genome_stratification_on_dataset_all:
    input:
        reference=BaseGenome.path(file_ending="/reference.fa.fai")
    output:
        stratification=GenomeStratificationOnDataset.path(stratification_type="all")
    run:
        with open(input.reference) as ref:
            with open(output.stratification, "w") as out:
                for line in ref:
                    l = line.split()
                    out.write(f"{l[0]}\t0\t{l[1]}\n")

ruleorder: subset_genome_stratification_on_dataset_all > subset_genome_stratification_on_dataset


rule subset_genotypes_on_stratification:
    input:
        vcf = GenotypeResults.path(),
        stratification=GenomeStratification.path()
    output:
        vcf = StratifiedGenotypeResults.path(file_ending="/genotypes.vcf"),
        gz = StratifiedGenotypeResults.path(file_ending="/genotypes.vcf.gz")
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools intersect -header -wa -f 0.5 -u -a {input.vcf} -b {input.stratification} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} &&"
        "bgzip -c {output.vcf} > {output.gz} "


# "all" is simply the vcf copied
rule all_stratification:
    input:
        vcf = GenotypeResults.path(),
    output:
        vcf = StratifiedGenotypeResults.path(stratification_type="all"),
        gz = StratifiedGenotypeResults.path(stratification_type="all", file_ending="/genotypes.vcf.gz")
    shell:
        "cat {input.vcf} | "
        "python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} && "
        "bgzip -c {output.vcf} > {output.gz}"


ruleorder: all_stratification > subset_genotypes_on_stratification


rule get_best_genotypes:
    input:
        vcf = StratifiedGenotypeResults.path(),
    output:
        vcf = BestGenotypes.path(),
        gz = BestGenotypes.path(file_ending="/genotypes.vcf.gz")
    shell:
        """
        python3 scripts/get_best_genotypes.py {input.vcf} {wildcards.ratio_of_best_genotypes} > {output.vcf} &&
        bgzip -c {output.vcf} > {output.gz}
        """



rule subset_individual_on_stratification:
    input:
        vcf = Individual.path(),
        stratification=GenomeStratification.path()
    output:
        vcf = StratifiedIndividual.path(),
        gz = StratifiedIndividual.path(file_ending="/individual.vcf.gz"),
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools intersect -header -wa -f 0.5 -u -a {input.vcf} -b {input.stratification} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} && "
        "bgzip -c {output.vcf} > {output.gz}"


# "all" is simply the vcf copied
rule all_stratification_individual:
    input:
        vcf = Individual.path(),
    output:
        vcf = StratifiedIndividual.path(stratification_type="all"),
        gz = StratifiedIndividual.path(stratification_type="all", file_ending="/individual.vcf.gz"),
    shell:
        "cat {input.vcf} "
        "| python3 scripts/filter_vcf_on_variant_type.py {wildcards.stratification_variant_type} > {output.vcf} && "
        "bgzip -c {output.vcf} > {output.gz}"



ruleorder: all_stratification_individual > subset_individual_on_stratification



rule filter_individual:
    input:
        individual = StratifiedIndividual.path(),
        population = PopulationWithoutIndividual.path(),
        population_index = PopulationWithoutIndividual.path(file_ending="/population_without_individual.vcf.gz.tbi")
    output:
        individual = FilteredIndividual.path(individual_filter="only_variants_in_population_isec"),
        gz = FilteredIndividual.path(individual_filter="only_variants_in_population_isec", file_ending="/individual.vcf.gz"),
    params:
        dir = lambda wildcards, input, output: "/".join(output.individual.split("/")[:-1])
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bgzip -c {input.individual} > {input.individual}.gz &&
        tabix -p vcf {input.individual}.gz &&
        bcftools isec -p {params.dir} {input.individual}.gz {input.population} -O z && 
        zcat {params.dir}/0002.vcf.gz > {output.individual} &&
        bgzip -c {output.individual} > {output.gz}
        """


rule filter_individual2:
    """
    Filter an individual so that it only keeps variants that are in the population
    """
    input:
        individual = StratifiedIndividual.path(),
        population = PopulationWithoutIndividual.path(),
    output:
        individual = FilteredIndividual.path(individual_filter="only_variants_in_population"),
        gz = FilteredIndividual.path(individual_filter="only_variants_in_population", file_ending="/individual.vcf.gz"),
        variant_ids = FilteredIndividual.path(individual_filter="only_variants_in_population", file_ending="/untypable_ids.tsv")
    shell:
        """
        python3 scripts/filter_vcf_not_in_other_vcf.py {input.population} {input.individual} {output.variant_ids} > {output.individual} &&
        bgzip -c {output.individual} > {output.gz}
        """



rule filter_individual_no_filter:
    input:
        individual = StratifiedIndividual.path(),
    output:
        individual = FilteredIndividual.path(individual_filter="none"),
        gz = FilteredIndividual.path(individual_filter="none", file_ending="/individual.vcf.gz"),
    shell:
        "cp {input.individual} {output.individual} && bgzip -c {output.individual} > {output.gz}"



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
        rm -rf {params.out_path} && 
        rtg format -o {params.out_path} {input}
        """


rule rtg_tools_vcf_eval:
    input:
        truth=FilteredIndividual.path(file_ending="/individual.vcf.gz"),
        truth_index=FilteredIndividual.path(file_ending="/individual.vcf.gz.tbi"),
        genotypes=StratifiedGenotypeResults.path(file_ending="/genotypes.vcf.gz"),
        genotypes_index=StratifiedGenotypeResults.path(file_ending="/genotypes.vcf.gz.tbi"),
        sdf=BaseGenome.path(file_ending="/sdf/done")
    output:
        report = VcfEvalReport.path(),
        roc = VcfEvalRoc.path()
    params:
        out_dir = lambda wildcards, input, output: "/".join(output.report.split("/")[:-1]) + "/rtg_tools_tmp",
        sdf = lambda wildcards, input, output: input.sdf.replace("/done", "")
    conda:
        "../envs/rtg_tools.yml"
    shell:
        """
        rm -rf {params.out_dir} &&
        
        rtg vcfeval -b {input.truth} -c {input.genotypes} -t {params.sdf} -o {params.out_dir} &&
        mv {params.out_dir}/summary.txt {output.report} &&
        #rtg rocplot {params.out_dir}/weighted_roc.tsv.gz --png {output.roc} && gio open {output.roc}
        rtg rocplot {params.out_dir}/weighted_roc.tsv.gz 
        """


rule happy_evaluation:
    input:
        truth=FilteredIndividual.path(file_ending="/individual.vcf.gz"),
        truth_index=FilteredIndividual.path(file_ending="/individual.vcf.gz.tbi"),
        genotypes=BestGenotypes.path(file_ending="/genotypes.vcf.gz"),
        genotypes_index=StratifiedGenotypeResults.path(file_ending="/genotypes.vcf.gz.tbi"),
        reference=BaseGenome.path()
    output:
        roc=VcfEvalRoc.path(file_ending="/happy/happy.roc.all.csv.gz")
    params:
        out_dir = lambda wildcards,input,output: "/".join(output.roc.split("/")[:-1]),
    conda:
        "../envs/happy.yml"
    shell:
        """
        mkdir -p {params.out_dir} &&
        echo {params.out_dir} && 
        hap.py {input.truth} {input.genotypes} -o {params.out_dir}/happy -r {input.reference} --roc GQ
        """


rule roc_plot:
    input:
        roc=VcfEvalRoc.path(file_ending="/happy/happy.roc.all.csv.gz")
    output:
        png=VcfEvalRoc.path(file_ending="/happy_roc.png")
    run:
        import gzip
        import pandas as pd
        import plotly.express as px
        import numpy as np

        with gzip.open(input[0], "rt") as f:
            f.readline()
            lines = f.readlines()
            lines = (l.split(",") for l in lines)
            lines = [line for line in lines if line[0] == "INDEL" and line[3] == "ALL" and line[6] != "*"]

            quality_scores = np.array([float(line[6]) for line in lines])
            recalls = np.array([float(line[7]) for line in lines])
            precisions = np.array([float(line[8]) for line in lines])

            sorting = np.argsort(-quality_scores)
            quality_scores = quality_scores[sorting]
            recalls = recalls[sorting]
            precisions = precisions[sorting]

            df = pd.DataFrame({"quality_score": quality_scores, "recall": recalls, "precision": precisions})
            print(df)
            fig = px.line(df, x="recall", y="precision", text="quality_score")

            # save plotly figure
            fig.write_image(output.png)

            #px.show()


rule pangenie_ids:
    """
    Get the ids.tsv file neede for running pangenie's evaluation script
    """
    input:
        vcf=FilteredIndividual.path(file_ending="/individual.vcf"),
        #untypable_ids="empty.tsv"  # don't need untypable, variants already filtered
    output:
        ids = FilteredIndividual.path(file_ending="/ids.tsv")
    shell:
        """
        #cat {input.vcf} | python3 scripts/pangenie_skip_untypeable.py {input} | 
        cat {input.vcf} | python3 scripts/pangenie_get_ids.py > {output}
        """


rule pangenie_genotype_accuracy:
    """
    Using pangenie's script
    """
    input:
        ids=FilteredIndividual.path(file_ending="/ids.tsv"),
        truth=FilteredIndividual.path(file_ending="/individual.vcf"),
        sample=BestGenotypes.path(file_ending="/genotypes.vcf"),
    output:
        report=touch(PangenieReport.path())
    shell:
        """
		python3 scripts/pangenie_genotype_evaluation.py {input.truth} {input.sample} {input.ids} --qual 0
        """




rule run_truvari:
    input:
        truth = FilteredIndividual.path(file_ending="/individual.vcf.gz"),
        truth_index = FilteredIndividual.path(file_ending="/individual.vcf.gz.tbi"),
        sample = BestGenotypes.path(file_ending="/genotypes.vcf.gz"),
        sample_index = BestGenotypes.path(file_ending="/genotypes.vcf.gz.tbi"),
        regions = GenomeStratificationOnDataset.path(),
        ref = BaseGenome.path()
    output:
        report = TruvariReport.path()
    params:
        out_dir = lambda wildcards, input, output: "/".join(output.report.split("/")[:-1])
    conda:
        "../envs/truvari.yml"
    shell:
        """
        rm -rf {params.out_dir} && 
        truvari bench \
        -b {input.truth} \
        -c {input.sample} \
        -o {params.out_dir} \
        --refdist 500 \
        --chunksize 500 \
        --reference {input.ref} \
        --pick ac \
        --includebed {input.regions} \
        --no-ref a && 
        mv {params.out_dir}/summary.json {output.report}
        """


rule process_truvari_report:
    input:
        report= TruvariReport.path()
    output:
        f1 = TruvariF1Score.path(),
        recall = TruvariRecall.path(),
        precision = TruvariPrecision.path(),
    run:
        import json
        with open(input.report) as f:
            report = json.load(f)
            f1_score = report["f1"]
            precision = report["precision"]
            recall = report["recall"]

            print(f"Precision: {precision}. Recall: {recall}. F1: {f1_score}")
            with open(output.f1, "w") as f:
                f.write(str(f1_score))

            with open(output.precision, "w") as f:
                f.write(str(precision))

            with open(output.recall, "w") as f:
                f.write(str(recall))