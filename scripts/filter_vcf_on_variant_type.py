import sys
import logging
logging.basicConfig(level=logging.INFO)
filter_type = sys.argv[1]


def get_type(ref_seq, alt_seq):
    alt_seqs = alt_seq.split(",")
    max_alt_seq_length = max([len(i) for i in alt_seqs])
    if len(ref_seq) == 1 and max_alt_seq_length == 1:
        return "snp"
    elif len(ref_seq) >= 50 or len(alt_seq) >= 50:
        return "sv"
    elif len(ref_seq) >= 10 or len(alt_seq) >= 10:
        return "large_indel"
    else:
        return "indel"


def keep_type(filter_type, variant_type):
    if filter_type == "all":
        return True
    elif filter_type == "svs" and variant_type == "sv":
        return True
    elif filter_type == "snps_indels" and variant_type in ("snp", "indel", "large_indel"):
        return True
    elif filter_type == "snps" and variant_type == "snp":
        return True
    elif filter_type == "indels" and variant_type in ("indel", "large_indel"):
        return True
    elif filter_type == "large_indels" and variant_type == "large_indel":
        return True
    return False


n_skipped = 0
n_total = 0
for line in sys.stdin:
    if line.startswith("#"):
        print(line.strip())
        continue

    variant_type = get_type(line.split()[3], line.split()[4])

    if keep_type(filter_type, variant_type):
        print(line.strip())
    else:
        n_skipped += 1
    n_total += 1

logging.info(f"Skipped {n_skipped}/{n_total} variants not matching type {filter_type}")
