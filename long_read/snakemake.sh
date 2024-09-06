#!/bin/bash

module load snakemake || exit 1

source myconda
mamba activate base
snakemake -pr --keep-going --retries 1 --latency-wait 120 --use-conda --use-envmodules -s workflow/snakefile.smk --profile workflow/snakemake_profile --rerun-triggers mtime 
