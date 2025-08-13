from macrodna import MaCroDNA
import pandas as pd
import gurobipy as gp
import numpy as np
import rpy2.robjects
from rpy2.robjects import pandas2ri
import multiprocessing as mp
from argparse import ArgumentParser, Namespace
import pickle
import random
from pathlib import Path

class MultiCellMatcher:
	outdir : Path
	random_seed : int
	gene_fraction : float
	infiles : dict[str,Path] = {}
	data : dict[str,pd.DataFrame] = {}
	barcode_transcribe_df : pd.DataFrame
	correlation_methods : list[str]
	samplesheet : Path
	gurobi_license_file : Path
	debug_mode : bool
	run_similarities : bool
	run_macrodna : bool

	def __init__(self,args:Namespace):
		self.outdir = Path(args.outdir)
		if not self.outdir.is_dir():
			self.outdir.mkdir(parents=True)
		self.random_seed = args.seed
		self.samplesheet = args.data_sheet
		self.correlation_methods = args.methods.split(",")
		self.gene_fraction = args.fraction
		self.barcode_transcribe_df = pd.read_csv(args.ground_truth,sep="\t",index_col="RNA")
		self.gurobi_license_file = args.license
		self.debug_mode = args.debug
		self.run_similarities = not args.skip_similarities
		self.run_macrodna = not args.skip_macrodna

		if not (self.run_similarities or self.run_macrodna):
			print("Please allow at least one computation")
			exit(1)

	def parse_samplesheet(self)->bool:
		try:
			with open(self.samplesheet,"r") as inhandle:
				#sample sheet structure: first line is header, all following lines are data
				# (currently only first three cols are used, which need to be data source and filepath)
				header = inhandle.readline()
				for line in inhandle:
					parts = line.strip().split(",")
					if len(parts) >= 2:
						data_src = parts[0].lower()
						data_file = Path(parts[1])
						self.infiles[data_src] = data_file
				
			if "dna" not in self.infiles:
				print("Please provide the DNA copynumber profle")
				return False
			if "rna" not in self.infiles:
				print("Please provide the raw gene expression counts")
				False
			return True
		except Exception as e:
			print(e)
			return False

	def prepare_data(self)->bool:
		try:
			#load all data
			rna_shared_cells = set(self.barcode_transcribe_df.index)
			for data_source in self.infiles:
				temp_df = pd.read_csv(self.infiles[data_source],sep="\t",index_col="ENSEMBL")
				self.data[data_source] = temp_df
				if data_source[0] != "dna":
					rna_shared_cells = rna_shared_cells.intersection(temp_df.columns)
			del temp_df
			#find shared cells
			shared_cells = pd.DataFrame({"RNA":list(rna_shared_cells)}).join(self.barcode_transcribe_df,on="RNA")
			shared_cells = shared_cells[shared_cells["DNA"].isin(self.data["dna"].columns)]
			print(f"Found {shared_cells.shape[0]} shared cells")

			#find genes shared between DNA and RNA
			shared_genes = pd.DataFrame(index=sorted(set(self.data["dna"].index).intersection(self.data["rna"].index)))
			total_genes = shared_genes.shape[0]
			#select random genes
			random.seed(self.random_seed)
			sel_genes = random.sample(list(range(total_genes)),int(total_genes*self.gene_fraction))
			shared_genes = shared_genes.iloc[sel_genes]
			print(f"Used {shared_genes.shape[0]} out of {total_genes} genes ({round(shared_genes.shape[0]/total_genes,4)*100}%)")
			if self.debug_mode:
				shared_genes.to_csv(self.outdir/f"replicate_{self.random_seed}_genes.tsv",sep="\t",index_label="ENSEMBL")
			print("Final shapes of dataframes are: ")
			for data_source in self.data:
				cell_col = "DNA" if data_source=="dna" else "RNA"
				self.data[data_source] = shared_genes.join(self.data[data_source],how="left").sort_index()[shared_cells[cell_col]]
				print(f"\t{data_source} : {self.data[data_source].shape}")
			return True
		except Exception as e:
			print(e)
			return False
	
	def find_matches(self,task:tuple[str,str])->pd.Series:
		compared_data_source,method = task
		#prepare matric for calculations
		merged_df = pd.concat([self.data["dna"],self.data[compared_data_source]],axis=1).astype(np.float32)
		#calulate similarities
		if method == "cosine":
			with pandas2ri.converter.context():
				rpy2.robjects.globalenv["profiles_df"] = merged_df.fillna(0)
			cosines = rpy2.robjects.r("lsa::cosine(as.matrix(profiles_df))")
			scores_df = pd.DataFrame(np.reshape(np.array(cosines),(len(merged_df.columns),len(merged_df.columns))),index=merged_df.columns,columns=merged_df.columns)
		elif method in ["pearson","spearman","kendall"]:
			scores_df = merged_df.corr(method=method,numeric_only=True)
		#save rwa scores for later analysis
		if self.debug_mode:
			scores_df.to_csv(self.outdir / f"replicate_{self.random_seed}_{compared_data_source}_{method}.tsv",sep="\t")
		#find best match for each
		scores_df = scores_df.loc[self.data["dna"].columns,self.data[compared_data_source].columns]
		best_score = scores_df.apply(lambda x:max(abs(x)),axis=0)

		is_best_match = scores_df.apply(lambda x : abs(x) == best_score,axis=1).T.replace(True,{col:col for col in scores_df.T.columns}).replace(False,pd.NA).T
		cell_matches = is_best_match.apply(lambda x:[k for k in x if not pd.isna(k)][0] if x.count() > 0 else "NA",axis=0)
		return pd.Series(cell_matches,name=f"{compared_data_source}_{method}").sort_index()

	def apply_macrodna(self)->bool:
		try:
			#construct gurobi env from file
			with open(self.gurobi_license_file,"r") as license_handle:
				gp_params = {}
				gp_params["WLSACCESSID"] = license_handle.readline().strip()
				gp_params["WLSSECRET"] = license_handle.readline().strip()
				gp_params["LICENSEID"] = int(license_handle.readline().strip())
				gp_env = gp.Env(params=gp_params)
			mdna = MaCroDNA(dna_df=self.data["dna"],rna_df=self.data["rna"],gp_env = gp_env)
			cell2cell,cell2cell_step = mdna.cell2cell_assignment()
			cell2cell_step.to_csv(self.outdir/f"replicate_{self.random_seed}_macrodna.tsv",sep="\t")
			with open(self.outdir/f"replicate_{self.random_seed}_macrodna.pkl","wb") as outhandle:
				pickle.dump(cell2cell_step,outhandle)
			return True
		except Exception as e:
			print(e)
			return False

	def calculcate_similarities(self)->bool:
		try:
			pool = mp.Pool()
			tasks = []
			for data_source in self.data:
				if data_source == "dna":
					continue
				for method in self.correlation_methods:
					tasks.append((data_source,method))
			print(f"Combinations to compute: {tasks}")
			assignments = pd.concat(list(pool.map_async(self.find_matches,tasks).get()),axis=1)
			assignments.to_csv(self.outdir/f"replicate_{self.random_seed}_similarity_assignments.tsv",sep="\t",index_label="RNA")
			return True
		except Exception as e:
			print(e)
			return False

	def match_cells(self)->None:
		success = self.parse_samplesheet()
		if not success:
			print("Problems parsing samplesheet")
			exit(1)
		success = self.prepare_data()
		if not success:
			print("Problems in data integration")
			exit(1)
		
		if self.run_macrodna:
			success = self.apply_macrodna()
			if not success:
				print("Problems with MaCroDNA")
				exit(1)
			print("Done with MaCroDNA")
		
		if self.run_similarities:
			success = self.calculcate_similarities()
			if not success:
				print("Problems with similarity calculations")
				exit(1)
			print("Done with similarity calculations")
		print("Done")

		
def main():
	p = ArgumentParser()
	p.add_argument("-g","--ground_truth",type=str,required=True,help="Path to groudn truth tsv")
	p.add_argument("-d","--data_sheet",type=str,required=True,help="Path to sample sheet containing the data sources and corresponding files")
	p.add_argument("-o","--outdir",type=str,required=False,default="./multi_matcher_outs",help="Path to output dir")
	p.add_argument("-s","--seed",type=int,required=False,default=1,help="Random seed to use for subsampling")
	p.add_argument("-f","--fraction",type=float,required=False,default=0.95,help="Fraction of genes to use")
	p.add_argument("-m","--methods",type=str,required=False,default="pearson,cosine,kendall,spearman",help="Comma-separated list of similarities to compute")
	p.add_argument("-l","--license",type=str,required=False,default="~/gurobi_license.txt",help="Path to a gurobi license file (only keys, one per line), only required when applying MaCroDNA")
	p.add_argument("--debug",action="store_true",help="Write out all pairwise scores for similarities, rather than just best matches")
	p.add_argument("--skip_similarities",action="store_true",help="Don't calculate similarities (equivalent to giving no methods)")
	p.add_argument("--skip_macrodna",action="store_true",help="Skip MaCroDNA step")
	args = p.parse_args()

	matcher = MultiCellMatcher(args)

	matcher.match_cells()
		
if __name__ == "__main__":
	main()