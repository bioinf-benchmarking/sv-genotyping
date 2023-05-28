import os

rule simulate_reads:
    input:
        reference=BaseGenome.path(),
        individual=Individual.path(),
    output:
        Reads.path()
    params:
        out_path = lambda wildcards, input, output: os.path.sep.join(output[0].split(os.path.sep)[:-1]),
    shell:
        """
        diploid_sequence_simulator \
        --coverage {wildcards.coverage} \
        --read-length {wildcards.read_length} \
        {input.individual} \
        {input.reference} \
        {params.out_path}
        """


rule uncompress_fastq:
    input:
        fastq=Reads.path()
    output:
        fastq=Reads.path(file_ending="/reads.fq")
    shell:
        """
        gunzip -c {input} > {output}
        """