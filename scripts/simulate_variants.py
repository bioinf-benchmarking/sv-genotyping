import logging
logging.basicConfig(level=logging.INFO)
from bionumpy.simulate.variants import simulate_variants
import bionumpy as bnp
import typer

app = typer.Typer()


@app.command()
def simulate(base_genome: str, snp_rate: float, small_indel_rate: float, sv_indel_rate: float, out_file_name: str):
    genome = bnp.Genome.from_file(base_genome)
    genome = genome.read_sequence(base_genome)

    with bnp.open(out_file_name + ".tmp.vcf", "w") as f:

        for variants in simulate_variants(genome,
                                          snp_rate,
                                          small_indel_rate,
                                          sv_indel_rate):
            # remove variants close to start of chrom
            # some genotypers (paragraph) don't like that
            variants = variants[variants.position >= 150]
            f.write(variants)

    # add header
    header = f"""##fileformat=VCFv4.1
##reference={base_genome}
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=TARGETPOS,Number=1,Type=String,Description="Target position for duplications.">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
"""
    for chromosome, size in genome.genome_context.chrom_sizes.items():
        header += f"##contig=<ID={chromosome},length={size}>\n"

    header += """#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT"""

    with open(out_file_name + ".tmp.vcf", "r") as f:
        with open(out_file_name, "w") as out_file:
            out_file.write(header)
            for line in f:
                out_file.write(line.strip() + "\tGT\n")


if __name__ == "__main__":
    app()
