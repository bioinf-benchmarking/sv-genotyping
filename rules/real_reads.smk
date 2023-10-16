

def get_real_reads_file(wildcards, input, output):
    sample_name = open(input.sample_name).read().strip()

    samples = open("thousand_genomes_samples.tsv")
    samples = [line.strip().split("\t") for line in samples]
    
    sample = [l for l in samples if l[5] == sample_name and l[3] == "alignment" and "cram" in l[0]]

    assert len(sample) == 1, f"Did not find sample {sample_name}"
    url = sample[0][0].replace("ftp://", "http://")
    print(url)

    return url



rule download_real_reads:
    input:
        sample_name = RawPopulation.path(file_ending="/unfiltered_population.random_sample_number_{individual_id}.txt")  # this file contains a random seeded sample name
    output:
        RealRawReads.path(file_ending="/real_reads.cram")
    params:
        file = get_real_reads_file
    shell:
        """
        wget {params.file} -O {output}
        """


rule get_real_reads_reference:
    output:
        "local_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell:
        "zcat local_data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz > {output}"


"""
These cram files need the specific grch38 reference used here to be decoded
"""
rule convert_cram_to_fq:
    input:
        reference="local_data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        reads=RealRawReads.path(file_ending="/real_reads.cram")
    output:
        RealRawReads.path(file_ending="/real_reads.fq.gz")
    threads:
        8
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools fastq --threads 8 --reference {input.reference} {input.reads} | gzip -c > {output}"

