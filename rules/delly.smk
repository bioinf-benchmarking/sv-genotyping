

rule run_delly_call:
    input:
        ref = BaseGenome.path(),
        #ref="human_g1k_v37_decoy.small.fasta",
        #samples=["mapped/a.bam"],
        alns = [Reads.path(file_ending="/reads.sorted.bam")],
        #index = [Reads.path(file_ending="/reads.sorted.bam.bai")]
        #index=["mapped/a.bam.bai"],
        #bed="test.bed.gz",  # optional
    output:
        vcf = GenotypeResults.path(method="delly",file_ending="/calls.vcf.gz"),
    threads: lambda wildcards: int(wildcards.n_threads)
    resources:
        mem_mb=4096,
    wrapper:
        "v3.0.0/bio/delly"


rule run_delly_genotype:
    input:
        vcf = GenotypeResults.path(method="delly",file_ending="/calls.vcf.gz"),
        ref = BaseGenome.path(),
        mapped_reads = Reads.path(file_ending="/reads.sorted.bam"),
    output:
        vcf= GenotypeResults.path(method="delly",file_ending="/genotypes.vcf.gz"),
    conda:
        "../envs/delly.yml"

    shell:
        "delly call -v {input.vcf} -o {output.vcf} -g {input.ref} {input.mapped_reads}"


rule decompress_delly:
    input:
        vcf = GenotypeResults.path(method="delly",file_ending="/genotypes.vcf.gz"),
    output:
        vcf= GenotypeResults.path(method="delly",file_ending="/genotypes.vcf"),
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools view {input.vcf} > {output.vcf}"