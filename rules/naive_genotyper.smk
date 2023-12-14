

rule run_naive_genotyper:
    input:
        variants = FilteredPopulationWithCollapsedAlleles.path()
    output:
        results=GenotypeResults.path(method="naive_genotyper"),
    benchmark:
        GenotypeResults.path(method="naive_genotyper", file_ending="/benchmark.csv")
    threads:
        lambda wildcards: int(wildcards.n_threads)
    shell:
        "kage naive_genotyper -p {input.variants} -o {output.results} "


rule run_naive_genotyper_with_glimpse:
    input:
        variants = FilteredPopulationWithCollapsedAlleles.path(),
        glimpse_chunks=directory(FilteredPopulationWithCollapsedAlleles.path(file_ending="/GLIMPSE_chunks")),
    output:
        results=GenotypeResults.path(method="naive_genotyper_with_glimpse"),
    benchmark:
        GenotypeResults.path(method="naive_genotyper_with_glimpse", file_ending="/benchmark.csv")
    threads:
        1
    shell:
        "kage naive_genotyper -p {input.variants} -o {output.results} --glimpse True --glimpse-chunks {input.glimpse_chunks} "
