#!/bin/bash
#SBATCH -J line2line_macrodna
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-X%1
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=
#SBATCH --mail-type=START,FAIL,END
#SBATCH -p 
#SBATCH -o logs/out/%x_%A_%a.out
#SBATCH -e logs/err/%x_%A_%a.err

# Fill in partition, email and replace X with the number of bootstrapping runs to create
# Modify outdir, profiles_dir, bootstrap_fraction and license as necessary

outdir = ../line2line_benchmarking
profiles_dir = ../line2line_profiles
bootstrap_fraction = 0.75
license = ~/license.txt

python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cn_shared.tsv \
 -r ${profiles_dir}/rna_cn_shared.tsv \
 -m macrodna\
 -o ${outdir}/rna \
 -f $bootstrap_fraction \
 -l $license \
 -s $SLURM_ARRAY_TASK_ID