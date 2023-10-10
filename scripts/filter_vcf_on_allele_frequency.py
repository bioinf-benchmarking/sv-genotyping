import gzip
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import bionumpy as bnp
import sys
import tqdm
import npstructures as nps

in_file = sys.argv[1]
sv_frequency = float(sys.argv[2])
snps_indels_frequency = float(sys.argv[3])

to_keep = []

with bnp.open(in_file, buffer_type=bnp.io.VCFBuffer) as f:
    for chunk in f.read_chunks(200000000):
        is_sv = (chunk.ref_seq.shape[1] >= 50) | (chunk.alt_seq.shape[1] >= 50)

        # hacky way to get allele frequency fast
        info_fields = bnp.io.strops.join(chunk.info, ";")
        info_fields = bnp.io.strops.split(info_fields, sep=";")
        af_strings = info_fields[bnp.str_equal(info_fields[:, 0:3], "AF=")]
        assert len(af_strings) == len(chunk), "Not all lines contain AF=?"

        #afs = np.array([float(af_string.to_string()) for af_string in af_strings[:, 3:]])
        afs = (af_string.to_string() for af_string in af_strings[:, 3:])
        afs = nps.RaggedArray(
            [
                [float(af) for af in af_string.split(",")]
                for af_string in afs
            ]
        )
        # use highest af
        afs = np.max(afs, axis=-1)

        keep = (is_sv & (afs >= sv_frequency)) | (~is_sv & (afs >= snps_indels_frequency))
        to_keep.append(keep)
        logging.info(f"Keeping {np.sum(keep)}/{len(keep)} variants. {np.sum(is_sv)} SVs, {np.sum(~is_sv)} SNPs/indels in chunk")

to_keep = np.concatenate(to_keep)

i = 0
for line in tqdm.tqdm(gzip.open(in_file, "rb"), total=len(to_keep)):
    line = line.decode("utf-8")
    if line.startswith("#"):
        print(line, end="")
        continue

    if to_keep[i]:
        print(line, end="")

    i += 1
