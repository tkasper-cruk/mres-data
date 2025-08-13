import pandas as pd
import scanpy as sc
import numpy as np
import rpy2.robjects
from rpy2.robjects import pandas2ri
import multiprocessing as mp
from argparse import ArgumentParser, Namespace
import pickle
import random
from pathlib import Path

def main():
	p = ArgumentParser()
	p.add_argument("-d","--dna",type=str,required=True,help="Path to DNA counts tsv")
	p.add_argument("-r","--rna",type=str,required=True,help="Path to RNA counts tsv")
	p.add_argument("-o","--outdir",type=str,required=False,default="./multi_matcher_outs",help="Path to output dir")
	p.add_argument("-s","--seed",type=int,required=False,default=1,help="Random seed to use for subsampling")
	p.add_argument("-f","--fraction",type=float,required=False,default=0.95,help="Fraction of cells to use")
	p.add_argument("-m","--methods",type=str,required=False,default="pearson,cosine,kendall,spearman",help="Comma-separated list of similarities to compute")
	p.add_argument("-l","--license",type=str,required=False,help="Path to gurobi license file")
	args = p.parse_args()

	outdir = Path(args.outdir)
	if not outdir.is_dir():
		outdir.mkdir(exist_ok=True,parents=True)
	dna_df = pd.read_csv(args.dna,sep="\t",index_col=0).sample(frac = args.fraction,axis=1,random_state=args.seed)
	rna_df = pd.read_csv(args.rna,sep="\t",index_col=0).sample(frac = args.fraction,axis=1,random_state=args.seed)
	print(dna_df.shape)
	print(rna_df.shape)
	merged_df = dna_df.join(rna_df,how="inner")
	print(merged_df.shape,merged_df.head(),sep='\n')

	for method in args.methods.split(","):
		if method == "cosine":
			with pandas2ri.converter.context():
				rpy2.robjects.globalenv["profiles_df"] = merged_df.fillna(0)
			cosines = rpy2.robjects.r("lsa::cosine(as.matrix(profiles_df))")
			scores_df = pd.DataFrame(np.reshape(np.array(cosines),(len(merged_df.columns),len(merged_df.columns))),index=merged_df.columns,columns=merged_df.columns)
		elif method in ["pearson","spearman","kendall"]:
			scores_df = merged_df.corr(method=method,numeric_only=True)
		elif method == "macrodna":
			if not args.license:
				print("Please provide a path to a gurobi license file")
				continue
			with open(args.license,"r") as license_handle:
				gp_params = {}
				gp_params["WLSACCESSID"] = license_handle.readline().strip()
				gp_params["WLSSECRET"] = license_handle.readline().strip()
				gp_params["LICENSEID"] = int(license_handle.readline().strip())
				gp_env = gp.Env(params=gp_params)
    		mdna = MaCroDNA(dna_df=dna_df,rna_df=rna_df,gp_env = gp_env)
    		_,scores_df = mdna.cell2cell_assignment()
			del _		
		else:
			continue
		
		if method != "macrodna":
			scores_df = scores_df.loc[dna_df.columns,rna_df.columns]
		scores_df.to_csv(Path(args.outdir) / f"{method}_{args.seed}.tsv",sep="\t")
		del scores_df
		print(method)
if __name__ == "__main__":
	main()