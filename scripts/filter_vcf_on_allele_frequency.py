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
    for i, chunk in enumerate(f.read_chunks(20000000)):
        logging.info("Processing chunk of size %d" % len(chunk))
        is_sv = (chunk.ref_seq.shape[1] >= 50) | (chunk.alt_seq.shape[1] >= 50)

        afs = chunk.info.AF[:, 0]
        #afs = bnp.io.strops.str_to_float(afs)

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
