# How to use this workflow
## 0 - Create the conda environment
* Create the conda environment used for all steps of this workflow from the provided .yaml file (`env.yaml`)
## 1 - Prepare the datasets
* Prepare the benchmarking dataset using the corresponding notebooks (in `notebooks/`)
* For cell to cell accuracy assessment, a tsv file containing the corresponding barcodes is required. The file for DEFND-seq/10X Epi Multiome is provided (`matched_barcodes.tsv`)
## 2 - Run cell assignment
* Use the provided slurm templates as guidance
* Since gurobi acadmic licenses only support few (1/2) parralel executions, it is recommended to split execution between MaCroDNA and all other processes
* The slurmscripts use relative paths, so they should be run from within `sbatch_scripts/`
## 3 - Assess accuracies
* Use the provided notebooks to assess the accuracies of the tested methods
* Since processing the score files takes a while and is highly parallelisable, it is recommended to use a jupyter server with a high number of available cores for this, since it speeds processing up dramatically