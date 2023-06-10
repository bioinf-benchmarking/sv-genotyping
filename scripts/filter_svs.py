import logging
logging.basicConfig(level=logging.INFO)
import gzip
import typer

app = typer.Typer()


@app.command()
def filter(vcf_file_name: str):

    n_kept = 0
    n_skipped = 0
    with gzip.open(vcf_file_name) as f:
        for line in f:
            line = line.decode("utf-8")
            if line.startswith("#"):
                print(line.strip())
                continue
            l = line.split()
            ref = l[3]
            alt = l[4]
            if "<" in ref or "<" in alt or ">" in ref or ">" in alt:
                n_skipped += 1
                continue
            print(line.strip())
            n_kept += 1

    logging.info(f"Skipped {n_skipped} SVs with unknown sequences")
    logging.info(f"Kept {n_kept} SVs with known sequences")


if __name__ == "__main__":
    app()