import logging
logging.basicConfig(level=logging.INFO)
import sys
from collections import defaultdict

if len(sys.argv) == 1 or sys.argv[1] == "-":
    vcf = sys.stdin
else:
    vcf = open(sys.argv[1])

n_snps = 0
n_indels_removed = 0

snp_positions = defaultdict(set)

# First get position of all SNPs
for i, line in enumerate(vcf):
    if line.startswith("#"):
        continue

    if i % 100000 == 0:
        logging.info("%d lines processed" % i)

    l = line.split()
    chromosome = l[0]
    position = int(l[1])

    if len(l[3]) == 1 and len(l[4]) == 1:
        n_snps += 1
        snp_positions[chromosome].add(position)

logging.info("Found %d snps" % n_snps)

if len(sys.argv) == 1 or sys.argv[1] == "-":
    vcf = sys.stdin
else:
    vcf = open(sys.argv[1])

indel_positions = defaultdict(set)

for i, line in enumerate(vcf):

    if i % 1000000 == 0:
        logging.info("%d lines processed" % i)

    if line.startswith("#"):
        print(line.strip())
        continue

    l = line.split()
    chromosome = l[0]
    pos = int(l[1])

    if "N" in l[3].upper() or "N" in l[4].upper():
        logging.info("Skipped variant %s,%s due to N in sequence" % (l[0], l[1]))
        continue

    if len(l[3]) != 1 or len(l[4]) != 1:  # indel
        size = len(l[3]) + 1

        overlapping = False
        for j in range(pos, pos + size):
            if j in snp_positions[chromosome]:
                overlapping = True

            # also check against indels included so far
            if j in indel_positions[chromosome]:
                overlapping = True

        for j in range(pos, pos + size):
            indel_positions[chromosome].add(j)

        if overlapping:
            n_indels_removed += 1
            continue

    print(line.strip())


logging.info("Indels removed: %d" % n_indels_removed)
