library(glue)
library(stringr)
library(data.table)
library(dplyr)
library(vcfR)
library(Matrix)
if (!require(numbat)){
  install.packages("numbat",dependencies=TRUE,repos="https://cloud.r-project.org/")
}
library(numbat)
system('export PATH="/tools/miniforge/envs/numbat_test/bin/:$PATH"')
system("echo $PATH")
args <- commandArgs(trailingOnly=TRUE)
mode <- '10X'
label <- args[1]
sample <- label
outdir <- paste0("/data/numbat_pileup/",label)
bam <- paste0("/data/cellranger/",label,"/outs/possorted_genome_bam.bam")
mtx_dir <- paste0("/data/cellranger/",label,"/outs/filtered_feature_bc_matrix")
barcode = paste0(mtx_dir,"/barcodes.tsv.gz")
ncores <- 2 * as.numeric(args[2])
gmap = "/ref/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"
eagle = "/ref/Eagle_v2.4.1/eagle"
cellsnp <- "/tools/miniforge/envs/numbat_test/bin/cellsnp-lite"
snpvcf = "/ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz"
paneldir = "/ref/1000G_hg19"
genome = "hg19"
message(paste0('Running in ', mode, ' mode'))
message(paste0('Using genome version: ', genome))
## check if files exist
if (!file.exists(bam)) {
  stop('BAM file not found')
}
if (!file.exists(snpvcf)) {
  stop('SNP VCF not found')
}
if (!file.exists(barcode)) {
  stop('Barcode file not found')
}
if (!file.exists(gmap)) {
  stop('Genetic map not found')
}
if (!file.exists(paneldir)) {
  stop('Phasing reference panel not found')
}
message(outdir,"\t",sample,"\t",glue('{outdir}/pileup/{sample}'))
dir.create(outdir, showWarnings = TRUE, recursive = TRUE)
dir.create(glue('{outdir}/pileup'), showWarnings = TRUE)
dir.create(glue('{outdir}/phasing'), showWarnings = TRUE)
dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = TRUE)
## pileup
cmd = glue(
  cellsnp,
  '-s {bam}',
  '-b {barcode}',
  '-O {outdir}/pileup/{sample}',
  '-R {snpvcf}',
  '-p {ncores}',
  '--minMAF 0',
  '--minCOUNT 2',
  '--UMItag Auto',
  '--cellTAG CB',
  .sep = ' ')
cmds = c(cmd)
cat('Running pileup\n')
script = glue('{outdir}/run_pileup.sh')
list(cmds) %>% fwrite(script, sep = '\n')
system(glue('chmod +x {script}'))
system(glue('sh {script} 2>&1 | tee {outdir}/pileup.log'), intern = FALSE)
## VCF creation
cat('Creating VCFs\n')
samples <- c(sample)
# read in the pileup VCF
vcfs = lapply(samples, function(sample) {
  vcf_file = glue('{outdir}/pileup/{sample}/cellSNP.base.vcf')
  if (file.exists(vcf_file)) {
    if (file.size(vcf_file) != 0) {
      vcf = vcfR::read.vcfR(vcf_file, verbose = F)
      if (nrow(vcf@fix) == 0) {
        stop(glue('Pileup VCF for sample {sample} has 0 variants'))
      }
      return(vcf)
    } else {
      stop('Pileup VCF is empty')
    }
  } else {
    stop('Pileup VCF not found')
  }
})
# Remove chr prefix if present
vcfs = lapply(vcfs, function(vcf){
  vcf@fix[,1] <- gsub("chr", "", vcf@fix[,1])
  return(vcf)
})
numbat:::genotype(label, samples, vcfs, glue('{outdir}/phasing'), chr_prefix = TRUE)
## phasing
cat('Running phasing\n')
eagle_cmd = function(chr) {
  paste(eagle,
        glue('--numThreads {ncores}'),
        glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'),
        glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'),
        glue('--geneticMapFile={gmap}'),
        glue('--outPrefix {outdir}/phasing/{label}_chr{chr}.phased'),
        sep = ' ')
}
cmds = lapply(1:22, function(chr){eagle_cmd(chr)})
script = glue('{outdir}/run_phasing.sh')
list(cmds) %>% fwrite(script, sep = '\n')
system(glue('chmod +x {script}'))
tryCatch({
  system2(script, stdout = glue("{outdir}/phasing.log"))
},
warning = function(w){
  stop('Phasing failed')
})
## Generate allele count dataframe
cat('Generating allele count dataframes\n')
if (genome == 'hg19') {
  gtf = gtf_hg19
} else {
  gtf = gtf_hg38
}
genetic_map = fread(gmap) %>%
  setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
  group_by(CHROM) %>%
  mutate(
    start = POS,
    end = c(POS[2:length(POS)], POS[length(POS)])
  ) %>%
  ungroup()
for (sample in samples) {
  # read in phased VCF
  vcf_phased = lapply(1:22, function(chr) {
    vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
    if (file.exists(vcf_file)) {
      fread(vcf_file, skip = '#CHROM') %>%
        rename(CHROM = `#CHROM`) %>%
        mutate(CHROM = str_remove(CHROM, 'chr'))
    } else {
      stop('Phased VCF not found')
    }
  }) %>%
    Reduce(rbind, .) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM)))
  pu_dir = glue('{outdir}/pileup/{sample}')
  # pileup VCF
  vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))
  # count matrices
  AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
  DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))
  cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)
  df = numbat:::preprocess_allele(
    sample = label,
    vcf_pu = vcf_pu,
    vcf_phased = vcf_phased,
    AD = AD,
    DP = DP,
    barcodes = cell_barcodes,
    gtf = gtf,
    gmap = genetic_map
  ) %>%
    filter(GT %in% c('1|0', '0|1'))
  fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
}
#cat('All done!\n')