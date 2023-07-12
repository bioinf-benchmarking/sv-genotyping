


rule genotype_accuracy:
    input:
        truth = Individual.path(),
        genotypes = GenotypeResults.path()
    output:
        recall = GenotypeRecall.path(),
        one_minus_precision = GenotypeOneMinusPrecision.path(),
        f1 = GenotypeF1Score.path(),
        report = GenotypeReport.path()
    script:
        "../scripts/genotype_accuracy.py"



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



