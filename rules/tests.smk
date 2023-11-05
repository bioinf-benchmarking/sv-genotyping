

rule test_kage:
    input:
        #f1 = "data/sacCer3/medium/0.01/0.001/0/0.3/0.8/25/1/150/10.0/kage/1/all/GenotypeF1Score.txt"
        #f1 = "data/sacCer3/medium/0.01/0.001/0/0.3/0.8/0.0001/0.0001/20/1/150/10.0/kage/1/all/GenotypeF1Score.txt"
        f1= "data/sacCer3/medium/0.01/0.001/0/0.3/0.8/from_pangenome/1/0.0001/0.0001/20/1.0/simulated/150/10.0/0.001/kage/1/all/all/none/1.0/all/GenotypeF1Score.txt"

    output:
        touch("test_kage")
    run:
        with open(input[0]) as f:
            f1 = float(f.read())
            assert f1 > 0.96

"""
        f1 = GenotypeF1Score.from_flat_params(
            genome_build="sacCer3",
            size="small",
            snp_rate = 0.05,
            small_indel_rate=0.01,
            sv_indel_rate = 0,
            n_individuals = 25,
            allele_frequency = 0.3,
            correlation = 0.8,
            individual_id="simulated_1",
            read_length=150,
            coverage=10.0,
            n_threads=1,
            method="kage").file_path()
        """



#rule test_kage_multiallelic: