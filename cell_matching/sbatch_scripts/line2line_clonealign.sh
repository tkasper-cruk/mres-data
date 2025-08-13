#!/bin/bash
#SBATCH -J line2line_clonealign
#SBATCH --cpus-per-task=8
#SBATCH --mem=450G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=
#SBATCH --mail-type=START,FAIL,END
#SBATCH -p 
#SBATCH -o /Users/kasper01/logs/out/%x_%J.out
#SBATCH -e /Users/kasper01/logs/err/%x_%J.err

# Fill in partition and email
# This requires large amounts of memory, 760 cells in cell-to-cell require ~ 450GB

outdir = ./line2line_benchmarking
profiles_dir = ./line2line_profiles

Rscript --vanilla ../scripts/clonealign.R \
 "${profiles_dir}/dna_cn_clonealign.tsv"\
 "${profiles_dir}/rna_cn_shared.tsv"\
 "${outdir}/clonealing_session.rds"\
 "${outdir}/clonealign_assignments.tsv"