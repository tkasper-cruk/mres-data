library(Seurat)
if (!require(copykat)){
    library(devtools)
    install_github("navinlabcode/copykat")
}
library(copykat)
library(biomaRt)

features_dir <- paste0(snakemake@input[["cr_dir"]],"/filtered_mtx")
message(features_dir)
message(getwd())
raw <- Read10X(data.dir=features_dir)
raw <- CreateSeuratObject(counts = raw, project = snakemake@params[["s_id"]], min.cells = 0, min.features = 0)
dir.create(snakemake@output[[1]])
setwd(snakemake@output[[1]])
exp.rawdata <- GetAssayData(raw,assay='RNA',layer = 'counts')
exp.rawdata <- exp.rawdata[order(rownames(exp.rawdata)),,drop=FALSE]
rn2 <- rownames(exp.rawdata)
mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl",host="http://grch37.ensembl.org"))
gene.ids <- getBM(filters = "ensembl_gene_id",attributes = c("ensembl_gene_id","hgnc_symbol"),values=rn2,mart = mart)
gene.ids <- gene.ids[order(gene.ids$ensembl_gene_id),,drop=FALSE]
gene.ids <- gene.ids[!duplicated(gene.ids$ensembl_gene_id),]
rownames(exp.rawdata) <- toupper(gene.ids$hgnc_symbol)
copykat.test <- copykat(
    rawmat=exp.rawdata, 
    id.type="S", 
    ngene.chr=5, 
    win.size=25, 
    KS.cut=0.1, 
    LOW.DR = 0.05,
    UP.DR = 0.2,
    sam.name=snakemake@params[["s_id"]], 
    distance=snakemake@config[["copykat_distance"]],
    norm.cell.names="",
    output.seg="FLASE",
    plot.genes="TRUE", 
    genome="hg20",
    n.cores=snakemake@threads)

save.image("copyKAT.RData")