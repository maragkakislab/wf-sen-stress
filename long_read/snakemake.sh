#!/bin/bash

module load snakemake || exit 1

source myconda
mamba activate base
snakemake -pr --keep-going --retries 5 --latency-wait 120 --use-conda --use-envmodules -s workflow/sra_snake.smk --profile workflow/snakemake_profile --rerun-triggers mtime 
