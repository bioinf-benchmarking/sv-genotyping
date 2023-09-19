import sys
import random

file = sys.argv[1]
seed = int(sys.argv[2])
n_lines = int(sys.argv[3])
out_file = sys.argv[4]

random.seed(seed)

with open(out_file, "w") as out:
    with open(file) as f:
        lines = list(f.readlines())
        random.shuffle(lines)

        out.writelines(
            [line.strip() + "\n" for line in lines][0:n_lines]
        )
