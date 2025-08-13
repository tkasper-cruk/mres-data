# How to use this workflow
## 1 - Set up the environment for first use
* Create the conda environment to compile the C++ scripts on your system (using the `resources/barcode_matching_source/compile_env.yaml` yaml file)
* Compile all C++ scripts for use using `resources/barcode_matching_source/compile_all.sh` in the created conda environment
* Modify workflow/profiles/sc_slurm/config.yaml to fit your environment and dataset size (the preset works well for ~1000 cells on two lanes of a 25B flow cell on NovaSeq X)
## 2 - Run the pipeline
* Modify config/main.yaml to fit your data
* Move the raw reads into `results/raw_fastq/{sample_id}`
* Run the pipeline using snakemake. Requires snakemake, the slurm executor plugin and the `--use-conda` flag. An example comand is in `run_workflow.sh`
* Using a screen session to avoid crashes due to client disconnects is recommended