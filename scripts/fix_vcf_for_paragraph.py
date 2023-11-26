"""
Adds a padding base to indels that do not start with the same base.
Paragraph requires this padding base, although it is not required by the VCF spec.
"""
import dataclasses
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import bionumpy as bnp

import sys
reference = sys.argv[1]
vcf = bnp.open(sys.argv[2])
out_file = bnp.open(sys.argv[3], "w")

genome = bnp.Genome.from_file(reference)
genome_sequence = genome.read_sequence()

for chunk in vcf.read_chunks():
    is_indel = (chunk.ref_seq.shape[1] > 1) | (chunk.alt_seq.shape[1] > 1)
    ref_bases = genome_sequence.extract_intervals(bnp.Interval(chunk.chromosome, chunk.position, chunk.position+1))
    ref_bases_before = genome_sequence.extract_intervals(bnp.Interval(chunk.chromosome, chunk.position-1, chunk.position))
    ref_bases = ref_bases[:, 0]
    ref_bases_before = ref_bases_before[:, 0]
    ref_mismatch = (ref_bases != chunk.ref_seq[:, 0]) | (ref_bases != chunk.alt_seq[:, 0])

    needs_padding = is_indel & ref_mismatch
    #print("Old")
    #print(chunk[needs_padding][0:10])

    ref_bases = ref_bases.to_string()
    ref_bases_before = ref_bases_before.to_string()

    new_ref = [
        seq.to_string() if not needs_padding[i] else ref_bases_before[i] + seq.to_string() for i, seq in enumerate(chunk.ref_seq)
    ]

    new_alt = [
        seq.to_string() if not needs_padding[i] else ref_bases_before[i] + seq.to_string() for i, seq in enumerate(chunk.alt_seq)
    ]

    logging.info(f"Padded {np.count_nonzero(needs_padding)}/{len(chunk)} variants. {np.sum(is_indel)} indels in chunk")

    new_position = chunk.position.copy()
    new_position[needs_padding] -= 1
    #assert np.all(chunk.position[~needs_padding] == new_position[~needs_padding])
    print(chunk.position[needs_padding])
    #print(chunk.position)
    #print(new_position)
    #print(chunk.position[needs_padding])
    #print(new_position[needs_padding])
    #print("..")

    chunk.position = new_position+1
    chunk.ref_seq = bnp.as_encoded_array(new_ref)
    chunk.alt_seq = bnp.as_encoded_array(new_alt)
    #print("New")
    #print(chunk[needs_padding][0:100])

    #chunk = dataclasses.replace(chunk, ref_seq=chunk.ref_seq, alt_seq=chunk.alt_seq, position=chunk.position)


    out_file.write(chunk)
