#!/bin/bash
#SBATCH -J cell2cell_clonealign
#SBATCH --cpus-per-task=8
#SBATCH --mem=450G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=
#SBATCH --mail-type=START,FAIL,END
#SBATCH -p 
#SBATCH -o logs/out/%x_%J.out
#SBATCH -e logs/err/%x_%J.err

# Fill in partition and email
# This requires large amounts of memory, 760 cells in cell-to-cell require ~ 450GB

outdir = ../cell2cell_benchmarking
profiles_dir = ../cell2cell_profiles

Rscript --vanilla clonealign.R \
 "${profiles_dir}/cell2cell_scabsolute_cn.tsv"\
 "${profiles_dir}/cell2cell_rna_cn.tsv"\
 "${outdir}/clonealing_session.rds"\
 "${outdir}/clonealign_assignments.tsv"