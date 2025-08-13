#!/bin/bash
#SBATCH -J line2line_similarities
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --array=1-X
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=
#SBATCH --mail-type=START,FAIL,END
#SBATCH -p 
#SBATCH -o logs/out/%x_%A_%a.out
#SBATCH -e logs/err/%x_%A_%a.err


# Fill in partition, email and replace X with the number of bootstrapping runs to create
# Modify outdir, profiles_dir and bootstrap_fraction as necessary

outdir = ./line2line_benchmarking
profiles_dir = ./line2line_profiles
bootstrap_fraction = 0.75

srun python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cn_allgenes.tsv \
 -r ${profiles_dir}/ckat_cn.tsv \
 -o ${outdir}/copykat \
 -f $bootstrap_fraction \
 -s $SLURM_ARRAY_TASK_ID &
srun python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cn_shared.tsv \
 -r ${profiles_dir}/rna_cn_shared.tsv \
 -o ${outdir}/rna \
 -f $bootstrap_fraction \
 -s $SLURM_ARRAY_TASK_ID &
srun python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cat_allgenes.tsv \
 -r ${profiles_dir}/numbat_cat.tsv \
 -o ${outdir}/numbat_cat \
 -f $bootstrap_fraction \
 -s $SLURM_ARRAY_TASK_ID &
srun python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cn_allgenes.tsv \
 -r ${profiles_dir}/numbat_cat.tsv \
 -o ${outdir}/numbat_cn \
 -f $bootstrap_fraction \
 -s $SLURM_ARRAY_TASK_ID &
srun python ../scripts/benchmark_line2line.py \
 -d ${profiles_dir}/dna_cn_allgenes.tsv \
 -r ${profiles_dir}/cvae_cn.tsv \
 -o ${outdir}/copyvae \
 -f $bootstrap_fraction \
 -s $SLURM_ARRAY_TASK_ID &

