import logging

import numpy as np

logging.basicConfig(level=logging.INFO)
import sys

"""
Filters a vcf by keeping a ratio of the best genotypes given QC scores
"""

ratio_to_keep = float(sys.argv[2])

assert ratio_to_keep >= 0.0
assert ratio_to_keep <= 1.0

with open(sys.argv[1]) as f:
    #lines = np.array([l for l in f.readlines() if not l.startswith("#")])
    #qualities = np.array([float(l[9].split(":")[-1]) for l in f if not l.startswith("#")])
    qualities = []
    for line in f:
        if not line.startswith("#"):
            try:
                qualities.append(float(line.split()[9].split(":")[-1]))
            except ValueError:
                logging.error("ERror parsing quality score from line. Setting score to 0")
                qualities.append(0)
                logging.error("Genotype column: %s" % line.split()[9])

    qualities = np.array(qualities)
    sorting = np.argsort(qualities)

    n_to_keep = int(len(qualities) * ratio_to_keep)
    split = len(qualities)-n_to_keep

    keep = set(sorting[split:])

    if ratio_to_keep == 1.0:
        assert len(keep) == len(qualities)

logging.info(f"Keeping {len(keep)} variants when picking the best {ratio_to_keep} of variants")


i = 0
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("#"):
            print(line.strip())
            continue

        if i in keep:
            print(line.strip())

        i += 1