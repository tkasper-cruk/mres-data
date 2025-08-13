#!/bin/bash
snakemake --use-conda  --rerun-incomplete --workflow-profile workflow/profiles/sc_slurm --profile slurm -p --slurm-keep-successful-logs
