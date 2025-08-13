# How to use this workflow
## 1 - Set up the environment
* Make sure to have cellranger and Eagle installed (set paths to the binary in the config file)
* Acquire the SNP reference and haplotype panels for numbat (see [https://kharchenkolab.github.io/numbat/articles/numbat.html](https://kharchenkolab.github.io/numbat/articles/numbat.html))
* Modify config/main.yaml to fit your data/requirements
* Modify workflow/profiles/sc_slurm/config.yaml to fit the data and cluster environment
* Move the raw reads into results/raw_reads/sample_id
## 2 - Run cellranger
* Run the workflow/Snakefile_cellranger file to run pre-processing, cellranger and cellbender (requires `--use-conda` flag)
* For use on a slurm cluster, use the run_cellranger.sh script (requires a snakemake user profile named slurm)
* Using a screen session is recommended to avoid crashes due to client disconnects
## 3 - Post-cellranger processing
### 3.1 Create QC conda env
* use the file in workflow/envs/cellranger_post.yaml to create a conda env containing all software required for the QC pipeline e.g. mamba env create -f workflow/envs/cellranger_post.yaml 
* this will create an environment called cellranger_qc, which should be used to runf the notebook
### 3.2 Run the notebook
* You will find the QC notebook in the results/cellranger/sample_id folder
* Run the notebook to run all QC/Filtering steps on the cellranger output
## 4 - Run scRNA-based CNV callers
* Run the workflow/Snakefile_cnv file to run pre-processing, cellranger and cellbender (requires `--use-conda` flag)
* For use on a slurm cluster, use the run_cnv.sh script (requires a snakemake user profile named slurm)
* using a screen session is recommended to avoid crashes due to client disconnects