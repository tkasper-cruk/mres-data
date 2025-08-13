#!/bin/bash
#SBATCH -J cell2cell_macrodna
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --array=1-X%1
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=
#SBATCH --mail-type=START,FAIL,END
#SBATCH -p 
#SBATCH -o logs/out/%x_%A_%a.out
#SBATCH -e logs/err/%x_%A_%a.err

# Fill in partition, email and replace X with the number of bootstrapping runs to create
# Modify outdir, profiles_dir, bootstrap_fraction, groundtruth and license as necessary

outdir = ../cell2cell_benchmarking
profiles_dir = ../cell2cell_profiles
bootstrap_fraction = 0.75
license = ~/license.txt
groundtruth = ../matched_barcodes.tsv

python ../scripts/assign_cells.py \
 -o $outdir \
 -g $groundtruth \
 -l $license \
 -d ${profiles_dir}/data_file.csv \
 -f $boostrap_fraction \
 --debug \
 --skip_similarities \
 -s $SLURM_ARRAY_TASK_ID