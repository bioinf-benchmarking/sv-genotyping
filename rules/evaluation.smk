


rule genotype_accuracy:
    input:
        truth = Individual.path(),
        genotypes = GenotypeResults.path()
    output:
        recall = GenotypeRecall.path(),
        one_minus_precision = GenotypeOneMinusPrecision.path(),
        f1 = GenotypeF1Score.path()
    script:
        "../scripts/genotype_accuracy.py"




