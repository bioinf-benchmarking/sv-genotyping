# import bionumpy as bnp
# genotypes = bnp.open(input.genotypes, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
# truth = bnp.open(input.truth, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
import logging
logging.basicConfig(level=logging.INFO)

def parse_genotype(genotype):
    if genotype == ".":
        return "0/0"
    else:
        return genotype.replace("|", "/").replace("1/0", "0/1")

def get_variant_type(ref, alt):
    if len(ref) == len(alt) == 1:
        return "snps"
    elif len(ref) <= 50 and len(alt) <= 50:
        return "small_indels"
    else:
        return "svs"


def get_genotypes_from_vcf(vcf_file_name, variant_type="all"):
    assert variant_type in ["all", "small_indels", "snps", "svs"]
    out = {}
    n_skipped = 0
    with open(vcf_file_name) as f:
        lines = (l for l in f if not l.startswith("#"))
        for line in lines:
            l = line.split()
            chrom = l[0]
            start = l[1]
            ref = l[3]
            alt = l[4]
            id = (chrom, start, ref, alt)
            genotype = l[9].split(":")[0]

            if variant_type != "all":
                if get_variant_type(ref, alt) != variant_type:
                    n_skipped += 1
                    continue

            out[id] = parse_genotype(genotype)
        #genotypes = (l.split()[9].split(":")[0] for l in lines)
        #return [parse_genotype(g) for g in genotypes]
        logging.info(f"Skipped {n_skipped} variants not matching variant type {variant_type}")
        return out


variant_type = snakemake.wildcards.variant_type

genotypes = get_genotypes_from_vcf(snakemake.input.genotypes, variant_type)
truth = get_genotypes_from_vcf(snakemake.input.truth, variant_type)

true_positive = 0
true_negative = 0
false_positive = 0
false_negative = 0

for id, t in truth.items():
    if id not in genotypes:
        # did not genotype, treat as 0/0
        g = "0/0"
    else:
        g = genotypes[id]

    # t = t.to_string().replace("|", "/").replace("1/0", "0/1")
    # g = g.to_string().replace("|", "/").replace("1/0", "0/1")

    if t == g and t != "0/0":
        true_positive += 1
    elif t == '0/0' and g != '0/0':
        false_positive += 1
    elif t != '0/0' and g != t:
        false_negative += 1
    elif t == "0/0" and g == "0/0":
        true_negative += 1
    else:
        assert False, (t, g)


for id, variant in genotypes.items():
    if id not in truth:
        logging.error("Found genotype not in truth set:")
        logging.error(f"{id}")
        raise Exception()

recall = true_positive / (true_positive + false_negative)
precision = true_positive / (true_positive + false_positive)
f1_score = 2 * (precision * recall) / (precision + recall)

print(f"Recall: {recall}, One minus precision: {1 - precision}, F1 score: {f1_score}")

with open(snakemake.output.recall, 'w') as f:
    f.write(str(recall))

with open(snakemake.output.one_minus_precision, 'w') as f:
    f.write(str(1 - precision))

with open(snakemake.output.f1, 'w') as f:
    f.write(str(f1_score))

