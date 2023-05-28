
def parse_genotype(genotype):
    if genotype == ".":
        return "0/0"
    else:
        return genotype.replace("|", "/").replace("1/0", "0/1")


def get_genotypes_from_vcf(vcf_file_name):
    out = {}
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
            out[id] = parse_genotype(genotype)
        #genotypes = (l.split()[9].split(":")[0] for l in lines)
        #return [parse_genotype(g) for g in genotypes]
        return out


rule genotype_accuracy:
    input:
        truth = Individual.path(),
        genotypes = GenotypeResults.path()
    output:
        recall = GenotypeRecall.path(),
        one_minus_precision = GenotypeOneMinusPrecision.path(),
        f1 = GenotypeF1Score.path()
    run:
        #import bionumpy as bnp
        #genotypes = bnp.open(input.genotypes, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
        #truth = bnp.open(input.truth, buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer).read()
        genotypes = get_genotypes_from_vcf(input.genotypes)
        truth = get_genotypes_from_vcf(input.truth)

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

            #t = t.to_string().replace("|", "/").replace("1/0", "0/1")
            #g = g.to_string().replace("|", "/").replace("1/0", "0/1")

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

        recall = true_positive / (true_positive + false_negative)
        precision = true_positive / (true_positive + false_positive)
        f1_score = 2 * (precision * recall) / (precision + recall)

        print(f"Recall: {recall}, One minus precision: {1 - precision}, F1 score: {f1_score}")

        with open(output.recall, 'w') as f:
            f.write(str(recall))

        with open(output.one_minus_precision, 'w') as f:
            f.write(str(1 - precision))

        with open(output.f1, 'w') as f:
            f.write(str(f1_score))



