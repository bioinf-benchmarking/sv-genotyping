import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import bionumpy as bnp
import typer
app = typer.Typer()

def simulate_population_genotype_matrix(n_variants: int, n_individuals: int, ratio_non_ref_allele: float, correlation: float):
    """
    Creates haplotypes.
    Allele frequencies are drawn from a geometric distribution with mean mean_allele_freq.
    """

    """
    ratio_non_ref_allele are ratio of variants where the most common allele is not the reference allele
    
    each haplotype follows the major allele with probability correlation_mean
    """

    n_haplotypes = n_individuals * 2
    major_allele = np.random.choice([0, 1], size=n_variants, p=[1 - ratio_non_ref_allele, ratio_non_ref_allele])
    follows_major = np.random.choice([0, 1], size=(n_haplotypes, n_variants), p=[1 - correlation, correlation])
    logging.info("N follows major: %s" % np.sum(follows_major, axis=1))
    haplotypes = (follows_major * major_allele).T

    logging.info(f"Haplotype matrix: {haplotypes}")
    logging.info("Allele frequencies: %s"  % (np.mean(haplotypes, axis=0)))

    # ensure at least one haplotype follows each allele
    allele_freq = np.sum(haplotypes, axis=1)
    followed_by_zero = (allele_freq == 0) | (allele_freq == n_haplotypes)
    logging.info(f"{np.sum(followed_by_zero)} variants are followed by all or none haplotypes. Swapping these at random haplotypes.")
    random_haplotypes = np.random.randint(0, n_haplotypes, size=np.sum(followed_by_zero))
    haplotypes[followed_by_zero==1, random_haplotypes] += 1
    haplotypes[followed_by_zero==1, random_haplotypes] %= 2

    genotype_matrix = np.zeros((n_variants, n_individuals))
    genotype_matrix += haplotypes[:, 0::2]*2
    genotype_matrix += haplotypes[:, 1::2]

    logging.info(f"Genotype matrix: {genotype_matrix}")
    logging.info("Allele frequencies genotype matrix: %s"  % (np.mean(haplotypes > 0, axis=1)))

    return genotype_matrix


@app.command()
def simulate_population_phased_vcf(vcf: str, out_file: str, n_individuals: int, ratio_non_ref_allele: float, correlation: float):
    """
    Creates a phased VCF file with n_individuals individuals.
    """

    variants = bnp.open(vcf).read()
    n_variants = len(variants)
    logging.info("Number of variants: %s" % n_variants)

    genotype_matrix = simulate_population_genotype_matrix(n_variants, n_individuals, ratio_non_ref_allele, correlation)
    genotype_strings = ["0|0", "0|1", "1|0", "1|1"]

    with open(vcf) as f:
        with open(out_file, "w") as out:
            variant_id = 0
            for line in f:
                if line.startswith("#"):
                    if line.lower().startswith("#chrom"):
                        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                        line = line.split()[0:9] + [f"simulated_{i}" for i in range(n_individuals)]
                        out.write("\t".join(line) + "\n")
                    else:
                        out.write(line)
                    continue

                fields = line.strip().split("\t")
                fields[8] = "GT"
                fields = fields[:9]
                for i in range(0, n_individuals):
                    fields.append(genotype_strings[int(genotype_matrix[variant_id, i])])

                out.write("\t".join(fields) + "\n")
                variant_id += 1


if __name__ == "__main__":
    app()
