# import bionumpy as bnp
# genotypes = bnp.open(input.genotypes, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
# truth = bnp.open(input.truth, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
import logging

import numpy as np

logging.basicConfig(level=logging.INFO)
from kage.analysis.genotype_accuracy import GenotypeAccuracy, IndexedGenotypes, IndexedGenotypes2, IndexedGenotypes3

config = snakemake.config["genomes"]

variant_type = snakemake.wildcards.limit_accuracy_to_variant_type

if len(snakemake.input) == 0:
    accuracy = np.nan
    recall = np.nan
    precision = np.nan
    f1_score = np.nan
    weighted_genotype_concordance = np.nan
    out_report = ""

else:

    input_is_biallelic = config[snakemake.wildcards.genome_build]["input_is_biallelic"]
    print("Input is biallelic? ", input_is_biallelic)

    if not input_is_biallelic:
        truth = IndexedGenotypes3.from_multiallelic_vcf(snakemake.input[0], convert_to_biallelic=False)
        sample = IndexedGenotypes3.from_multiallelic_vcf(snakemake.input[1], convert_to_biallelic=False)
        sample.normalize_against_reference_variants(truth)
    else:
        truth = IndexedGenotypes2.from_biallelic_vcf(snakemake.input[0])
        sample = IndexedGenotypes2.from_biallelic_vcf(snakemake.input[1])

    accuracy = GenotypeAccuracy(truth, sample, limit_to=variant_type)
    recall = accuracy.recall()
    precision = accuracy.precision()
    f1_score = accuracy.f1()
    weighted_genotype_concordance = accuracy.weighted_concordance
    weighted_genotype_concordance_pangenie_definition = accuracy.weighted_concordance_pangenie_definition
    out_report = accuracy.get_debug_report()

    print(f"Recall: {recall}, One minus precision: {1 - precision}, F1 score: {f1_score}, Weighted concordance: {weighted_genotype_concordance}. Weighted concordance pangenie definition: {weighted_genotype_concordance_pangenie_definition}")

with open(snakemake.output.recall, 'w') as f:
    f.write(str(recall))

with open(snakemake.output.one_minus_precision, 'w') as f:
    f.write(str(1 - precision))

with open(snakemake.output.f1, 'w') as f:
    f.write(str(f1_score))

with open(snakemake.output.report, 'wb') as f:
    pickle.dump(out_report, f)

