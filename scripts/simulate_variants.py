from bionumpy.simulate.variants import simulate_variants
import bionumpy as bnp

genome = bnp.Genome.from_file(snakemake.input.base_genome)
genome = genome.read_sequence(snakemake.input.base_genome)

with bnp.open(snakemake.output.variants + ".tmp.vcf", "w") as f:

    for variants in simulate_variants(genome,
                                      float(snakemake.wildcards.snp_rate),
                                      float(snakemake.wildcards.small_indel_rate),
                                      float(snakemake.wildcards.sv_indel_rate)):
        print(variants)
        f.write(variants)
        
        
# add header

header = f"""##fileformat=VCFv4.1
##reference={snakemake.input.base_genome}
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

with open(snakemake.output.variants + ".tmp.vcf", "r") as f:
    with open(snakemake.output.variants, "w") as out_file:
        out_file.write(header)
        for line in f:
            out_file.write(line.strip() + "\tGT\n")
