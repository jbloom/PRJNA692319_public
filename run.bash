#!/bin/bash
#SBATCH --cpus-per-task=16

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

printf "Running snakemake...\n"
snakemake -j 16 --keep-going --restart-times 3 --latency-wait 60 --rerun-incomplete --use-conda
printf "Run of snakemake complete.\n"
