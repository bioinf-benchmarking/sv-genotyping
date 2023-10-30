import logging
import bionumpy as bnp
import sys
from isal import igzip

import numpy as np
from tqdm import tqdm

base = sys.argv[1]
other = sys.argv[2]
out_ids = sys.argv[3]  # writes ids of removed variants to this file


logging.basicConfig(level=logging.INFO)
keep = set()


def get_id_from_line(line):
    l = line.split("\t")
    return l[0], l[1], l[3], l[4]


with igzip.open(base, "rb") as f:
    for line in tqdm(f):
        line = line.decode("utf-8")
        if line.startswith("#"):
            continue
        id = get_id_from_line(line)
        assert id not in keep, "Duplicate variant"
        keep.add(id)

n_skipped = 0
n_tot = 0
with open(out_ids, "w") as out_ids:
    with open(other) as f:
        for line in tqdm(f):
            if line.startswith("#"):
                print(line.strip())
                continue

            if get_id_from_line(line) in keep:
                print(line.strip())
            else:
                n_skipped += 1
                # pangenies code:
                fields = line.split()
                info_fields = {i.split('=')[0]: i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
                assert 'ID' in info_fields
                var_ids = info_fields['ID'].split(',')
                var_id = var_ids[0]

                out_ids.write(var_id + "\n")
                #logging.info("Skipping %s" % str(get_id_from_line(line)))

            n_tot += 1

logging.info(f"Skipped {n_skipped} of {n_tot} variants")


"""
all_filters = np.concatenate(all_filters)

with bnp.open(other) as f:
    i = 0
    for line in tqdm(f, total=len(all_filters)):
        if line.startswith("#"):
            print(line.strip())
            continue

        if all_filters[i]:
            print(line.strip())
        i += 1
"""