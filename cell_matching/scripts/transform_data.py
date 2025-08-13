import scanpy as sc
import numpy as np
import pandas as pd
from argparse import ArgumentParser, Namespace
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from sys import argv
from pathlib import Path

def apply_cat_mask(cn:int)->int:
    if cn < 2:
        return 1
    if cn > 2:
        return 3
    return 2

def main(args:list[str])->None:
	tools = ["copykat","copyvae","cellranger","cellbender","scabsolute","scanpy"]

	p = ArgumentParser()
	p.add_argument("-i","--input",type=str,required=True,help="Path to input direcory")
	p.add_argument("-n","--name",type=str,required=True,help="Sample name")
	p.add_argument("-o","--output",type=str,required=True,help="Path to output directory")
	p.add_argument("-t","--tool",type=str,required=True,help="Tool type")
	p.add_argument("-g","--gene_reference",type=str,required=False,default="",help="Path to gene annotations, only for scAbsolute or numbat")
	p.add_argument("-m","--mode",type=str,required=False,default="cn",help="Produce integer copy number (cn,default) or state (cat) profiles")

	args = p.parse_args(args)
	outdir = Path(args.output)
	if not outdir.is_dir():
		outdir.mkdir(parents=True)
	t = args.tool.lower()
	if args.tool.lower() not in tools:
		print(f"Tool {args.tool} not recognised, please choose from the following list: {' '.join([tools])}")
		exit(1)

	if t == "copyvae":
		#transform copyvae output
		copyvae = np.load(f"{args.input}/copy.npy")
		cvae_adata = sc.read_h5ad(f"{args.input}/filtered.h5ad")
		transformed_df = pd.DataFrame(index=cvae_adata.var["gene_ids"],columns=cvae_adata.obs_names).rename(columns=lambda x:x[:-2])
		for gene_id in range(copyvae.shape[1]):
			for cell_id in range(copyvae.shape[0]):
				transformed_df.iloc[gene_id*25:gene_id*25+25,cell_id] = copyvae[cell_id,gene_id]
	
	elif t == "copykat":
		infile = Path(args.input) / f"{args.name}_copykat_CNA_raw_results_gene_by_cell.txt"
		copykat_df = pd.read_csv(infile,sep="\t")
		copykat_df = copykat_df.set_index("ensembl_gene_id").drop(columns=["abspos","chromosome_name","start_position","end_position","hgnc_symbol","band"])
		transformed_df = copykat_df.rename(columns=lambda x:x[:-2])
	
	elif t in ["scanpy","rna"]:
		infile = Path(args.input) / "filtered.h5ad"
		rna_data = sc.read_h5ad(infile)
		rna_data.var_names = rna_data.var["gene_ids"]
		transformed_df = rna_data.to_df().T.rename(columns=lambda x: x[:-2])

	elif t == "numbat":
		if args.mode != "cat":
			print("Extracting integer profiles for numbat is not possible")
			exit(1)
		robjects.globalenv["indir"] = args.input
		robjects.globalenv["genefile"] = args.gene_reference
		robjects.r(
			"""
			library(tidyverse)
			replace_cat <- function(x){
				if (x == "del"){
					return(1)
				}
				if (x == "neu"){
					return(2)
				}
				if (x == "amp"){
					return(3)
				}
				return(NA)
				}
			cells <- read_tsv(paste0(indir,"/clone_post_2.tsv")) %>%
			select(c(cell,clone_opt)) %>%
			mutate(
				cell = str_split_i(cell,"-",1)
			)
			segs <- read_tsv(paste0(indir,"/bulk_clones_final.tsv.gz")) %>%
			distinct(CHROM,seg,seg_start,seg_end,cnv_state,sample) %>%
			mutate(
				cnv_state = as.numeric(map(cnv_state,replace_cat))
			)
			genes_per_cell <- read_tsv(
				genefile,
				col_names = c("CHROM","START","END","HGNC","ENSEMBL","arm","loc","band")) %>%
				filter(
					CHROM != "X",
					CHROM != "Y"
				) %>%
			mutate(
				CHROM = as.numeric(CHROM)
			) %>%
			select(c(1,2,3,5)) %>%
			left_join(
				segs,
				by = join_by(
				CHROM == CHROM,
				within(START,END,seg_start,seg_end)
				)
			) %>%
			drop_na()  %>%
			select(c(ENSEMBL,cnv_state,sample)) %>%
			left_join(
				cells,
				by=join_by(sample==clone_opt),
				relationship = "many-to-many"
				) %>%
			select(c(-sample)) %>%
			pivot_wider(
				names_from = cell,
				values_from = cnv_state
			)
			"""
		)
		with pandas2ri.converter.context():
			transformed_df = robjects.globalenv["genes_per_cell"]
		transformed_df = transformed_df.set_index("ENSEMBL").dropna().rename(columns=lambda x: x.replace("_","-"))

	elif t == "scabsolute":
		if not args.gene_reference:
			print("Please provide gene reference")
			exit(1)
		robjects.globalenv["infile"] = args.input
		robjects.globalenv["genefile"] = args.gene_reference
		robjects.r(
			"""
			library(tidyverse)
			object <- readRDS(infile)
			dna_cn <- object@assayData$copynumber
			bins <- rownames(dna_cn)
			dna_cn <- as_tibble(dna_cn) %>% 
			mutate(
				CHR = str_split_i(bins,":",1),
				START = as.numeric(str_split_i(str_split_i(bins,":",2),"-",1)),
				END = as.numeric(str_split_i(bins,"-",2))
			)

			gene_cn <- read_delim(
				file = genefile,
				delim = "\t",
				col_names = c("CHR","START","END","HGNC","ENSEMBL","ARM","LOC","BAND")
			) %>%
			full_join(
				dna_cn,
				by = join_by(CHR, overlaps(START,END,START,END))
			) %>% 
			select(c(-START.x,-START.y,-END.x,-END.y,-ARM,-BAND,-LOC,-CHR,-HGNC))
			genes <- list(gene_cn$ENSEMBL)
			gene_cn <- as_tibble(
			aggregate(as.data.frame(select(gene_cn,-1)),by=genes,FUN = median)) %>%
			rename(ENSEMBL=Group.1)
			"""
			)
		with pandas2ri.converter.context():
			transformed_df = robjects.globalenv["gene_cn"]
		transformed_df = transformed_df.set_index("ENSEMBL").dropna().rename(columns=lambda x: x.replace("_","-"))
		if args.mode == "cat":
			transformed_df = transformed_df.map(apply_cat_mask)
		

	transformed_df.to_csv(outidr/f"{args.name}_{args.tool}_{args.mode}.tsv",sep="\t",index_label="ENSEMBL")

if __name__ == "__main__":
	main(argv[1:])
