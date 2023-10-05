

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
            Individual.path(file_ending="/individual.vcf"),
            GenotypeResults.path()
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
        report = GenotypeReport.path()
    script:
        "../scripts/genotype_accuracy.py"



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



