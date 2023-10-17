import os

# GLIMPSE has no working conda or anything
# to make installation backwards compatible, download static binaries from release
rule install_glimpse:
    output:
        "glimpse/GLIMPSE_chunk_static",
        "glimpse/GLIMPSE_phase_static",
        "glimpse/GLIMPSE_ligate_static",
    shell:
        "wget -O glimpse/GLIMPSE_chunk_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static && "
        "wget -O glimpse/GLIMPSE_phase_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static && "
        "wget -O glimpse/GLIMPSE_ligate_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_ligate_static && "
        "chmod a+x glimpse/GLIMPSE_*_static "


rule make_glimpse_chunks:
    input:
        glimpse_command="glimpse/GLIMPSE_chunk_static",
        variants = FilteredPopulation.path(),
        variants_index = FilteredPopulation.path(file_ending="/filtered_population.vcf.gz.tbi"),
    output:
        FilteredPopulation.path(file_ending="/glimpse_chunks.txt"),
        #"data/{dataset}/glimpse_chunks.txt"
    params:
        regions = lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"],
        base_path = lambda wildcards, input, output: os.sep.join(output[0].split(os.sep)[:-1]),
    shell:
        """
        chromosomes='{params.regions}'
        for chromosome in $(echo $chromosomes | tr "," "\n")
            do
            {input.glimpse_command} --input {input.variants} --region $chromosome  --window-size 2000000 --buffer-size 200000 --output {params.base_path}/glimpse_chunk.$chromosome.txt &
        done
        wait
        cat {params.base_path}/glimpse_chunk.*.txt > {output}
        """


rule run_glimpse:
    input:
        vcf = GenotypeResults.path(method="kage_no_imputation"),
        ref_vcf = FilteredPopulation.path(),
        #vcf="data/{dataset}/kageNoPriorsN{n_individuals}all_{experiment}.vcf.gz",
        #ref_vcf="data/{dataset}/variants_{n_individuals}all.vcf.gz",
        glimpse_command="glimpse/GLIMPSE_phase_static",
        glimpse_command_ligate="glimpse/GLIMPSE_ligate_static",
        chunks=FilteredPopulation.path(file_ending="/glimpse_chunks.txt"),
        #chunks="data/{dataset}/glimpse_chunks.txt"
    output:
        vcf = GenotypeResults.path(method="kage_with_glimpse"),
        #vcf="data/{dataset}/kageWithGlimpseN{n_individuals,\d+}all_{experiment}.vcf.gz"
    benchmark:
        GenotypeResults.path(method="kage_with_glimpse", file_ending="/glimpse.csv")
        #"data/{dataset}/benchmarks/glimpse_{n_individuals,\d+}all_{experiment}.tsv"
    conda: "../envs/bcftools.yml"
    threads: lambda wildcards: int(wildcards.n_threads)
    params:
        regions = lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.size]["chromosomes"],
        base_path = lambda wildcards, input, output: os.sep.join(output[0].split(os.sep)[:-1]),
    shell:
        """
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf -f {input.vcf}.gz
        rm -f {params.base_path}/GLIMPSE-*.bcf
        cat {input.chunks} | parallel -j {config[n_threads]} --line-buffer "scripts/run_glimpse.sh {{}} {input.vcf}.gz {input.ref_vcf}"
        mkdir -p {params.base_path}

        # merge all result files 
        rm -rf {params.base_path}/glimpse_tmp_*.vcf.gz
        chromosomes='{params.regions}'
        for chromosome in $(echo $chromosomes | tr "," "\n")
            do
            LST={params.base_path}/glimpse_list$chromosome.tmp.txt
            ls {params.base_path}
            ls {input.vcf}.gz-GLIMPSE-$chromosome.*.bcf | python3 scripts/sort_glimpse_list.py > $LST
            {input.glimpse_command_ligate} --input $LST --output {params.base_path}/glimpse_tmp_$chromosome.vcf.gz
            tabix -p vcf -f {params.base_path}/glimpse_tmp_$chromosome.vcf.gz
        done
        wait
        bcftools concat {params.base_path}/glimpse_tmp_*.vcf.gz > {output.vcf}
        #tabix -p vcf -f {output.vcf}

        """
