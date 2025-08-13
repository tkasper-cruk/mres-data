tensorflow::use_condaenv("merge_dna_rna")
if (!require(clonealign)){
  devtools::install_github("kieranrcampbell/clonealign")
  library(clonealign)
}
args = commandArgs(trailingOnly=TRUE)


dna_data <- read.csv(args[1],sep="\t")
rna_data <- read.csv(args[2],sep="\t")
rownames(rna_data) <- rna_data$ENSEMBL
rna_data$ENSEMBL <- NULL
rna_mat <- as.matrix(rna_data)
rownames(dna_data) <- dna_data$ENSEMBL
dna_data$ENSEMBL <- NULL
dna_data <- dna_data + 0.01

cal <- clonealign(t(rna_mat),dna_data)
saveRDS(cal,args[3])

#cal <- readRDS("cal_object.RDS")

cal_assignments <- data.frame(cal$clone)
cal_assignments$RNA <- substr(colnames(rna_data),1,16)
cal_assignments$DNA <- substr(cal_assignments$cal.clone,stringr::str_length(cal_assignments$cal.clone)-15,stringr::str_length(cal_assignments$cal.clone)+1)
cal_assignments$cal.clone <- NULL
write.csv(cal_assignments,args[4],sep="\t")