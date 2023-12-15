
### Reproducing the KAGE2 experiments

Start by following the guide in the Readme of this repository for setting up the benchmarking pipeline.

Reproducing the figures presented in the KAGE2 paper can then be done by running (one command for each figure):

```bash
snakemake plots/figure1a.png --use-conda
snakemake plots/figure1b.png --use-conda
snakemake plots/figure2.png --use-conda
snakemake plots/figure3a.png --use-conda
snakemake plots/figure3b.png --use-conda
snakemake plots/figure4.png --use-conda

snakemake plots/supplementary_figure1.png --use-conda
snakemake plots/supplementary_figure2.png --use-conda

```