import os

rule simulate_reads:
    input:
        reference=BaseGenome.path(),
        index=BaseGenome.path(file_ending="/reference.fa.fai"),
        individual=Individual.path(),
    output:
        Reads.path(file_ending="/reads.fq.gz")
    threads:
        4
    params:
        out_path = lambda wildcards, input, output: os.path.sep.join(output[0].split(os.path.sep)[:-1]),
    shell:
        """
        diploid_sequence_simulator \
        --coverage {wildcards.coverage} \
        --read-length {wildcards.read_length} \
        --random-seed {config[random_seed]} \
        --n-threads 4 \
        --file-ending .fq.gz \
        {input.individual} \
        {input.reference} \
        {params.out_path}
        """


"""
rule compress_fastq:
    input:
        Reads.path(file_ending="/reads.fq")
    output:
        Reads.path(file_ending="/reads.fq.gz")
    shell:
        "gzip -c {input} > {output}"
"""

rule uncompress_fastq:
    input:
        fastq=Reads.path()
    output:
        fastq=Reads.path(file_ending="/reads.fq")
    shell:
        "gunzip -c {input} > {output}"