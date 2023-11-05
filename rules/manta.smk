

rule manta:
    input:
        ref = BaseGenome.path(),
        #ref="human_g1k_v37_decoy.small.fasta",
        #samples=["mapped/a.bam"],
        samples = [Reads.path(file_ending="/reads.sorted.bam")],
        index = [Reads.path(file_ending="/reads.sorted.bam.bai")]
        #index=["mapped/a.bam.bai"],
        #bed="test.bed.gz",  # optional
    output:
        vcf = GenotypeResults.path(method="manta", file_ending="/genotypes.vcf.gz"),
        idx = GenotypeResults.path(method="manta", file_ending="/genotypes.vcf.gz.tbi"),
        #vcf="results/out.bcf",
        #idx="results/out.bcf.csi",
        #cand_indel_vcf="results/small_indels.vcf.gz",
        cand_indel_vcf = GenotypeResults.path(method="manta", file_ending="/small_indels.vcf.gz"),
        #cand_indel_idx="results/small_indels.vcf.gz.tbi",
        cand_indel_idx = GenotypeResults.path(method="manta", file_ending="/small_indels.vcf.gz.tbi"),
        #cand_sv_vcf="results/cand_sv.vcf.gz",
        cand_sv_vcf = GenotypeResults.path(method="manta", file_ending="/cand_sv.vcf.gz"),
        #cand_sv_idx="results/cand_sv.vcf.gz.tbi",
        cand_sv_idx = GenotypeResults.path(method="manta", file_ending="/cand_sv.vcf.gz.tbi"),
    params:
        extra_cfg="",  # optional
        extra_run="",  # optional
    threads: lambda wildcards: int(wildcards.n_threads)
    resources:
        mem_mb=4096,
    wrapper:
        "v2.9.0/bio/manta"