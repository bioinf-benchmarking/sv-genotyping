
rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("{genome}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.2/bio/bwa/index"


rule bwa_map:
    input:
        reads = Reads.path(),
        idx = BaseGenome.path(file_ending=["/reference.fa", "/reference.fa.amb","/reference.fa.ann","/reference.fa.bwt","/reference.fa.pac","/reference.fa.sa"])
    output:
        reads=Reads.path(file_ending="/reads.bam"),
        sorted=Reads.path(file_ending="/reads.sorted.bam"),
        sorted_idx=Reads.path(file_ending="/reads.sorted.bam.bai")
    threads: config["n_threads"]
    conda: "../envs/bwa.yml"
    shell:
        """
        bwa mem -p -t {config[n_threads]} -R "@RG\\tID:sample\\tSM:sample" {input.idx[0]} {input.reads} | samtools view -o {output.reads} - && \
        sambamba sort {output.reads} && sambamba index -p {output.sorted}
        """

