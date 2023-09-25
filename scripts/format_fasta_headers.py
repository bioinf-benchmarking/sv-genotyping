"""
Removes coordinate from fasta headers
"""
import bionumpy as bnp
import sys
from dataclasses import replace

in_file = sys.argv[1]
out_file = sys.argv[2]

with bnp.open(in_file) as f:
    data = f.read()
    names = data.name
    new_names = bnp.as_encoded_array([name.to_string().split(":")[0] for name in names])
    data = replace(data, name=new_names)

    with bnp.open(out_file, "w") as out:
        out.write(data)
