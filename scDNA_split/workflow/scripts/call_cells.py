#general requirements
from argparse import ArgumentParser
from pathlib import Path
#data handling
import polars as pl
#otsu filter
from skimage.filters import threshold_otsu
#plotting
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def get_parser()->ArgumentParser:
    p = ArgumentParser()
    p.add_argument("-b","--barcode_filestem",required=True,help="The barcode coutns files to use as a basis for cell calling")
    p.add_argument("-o","--output_dir",type=str,required=True,help="Directory to store the selected barcodes in")
    p.add_argument("-s","--stats_dir",type=str,required=True,help="Directory to store the full table and plots in")
    p.add_argument("-n","--num_files",type=int,required=True,help="The number of ID files to process (numbered with -p digits, 0-based counting)")
    p.add_argument("-p","--padding",type=int,required=False,default=6,help="The number of digits to pad numbers to")
    p.add_argument("--slx",type=str,required=True,help="The SLX-ID of the run")
    p.add_argument("--sample",type=str,required=True,help="The sample ID")
    p.add_argument("--flowcell",type=str,required=True,help="The Flowcell ID")
    p.add_argument("--lane",type=str,required=True,help="The lane number")
    p.add_argument("--min_readcount",type=int,required=True,help="The minimum number of reads per cell to require")
    p.add_argument("--ordmag_expected",type=int,required=False,default=1,help="Number of expected cells for OrdMag thresholding")
    return p

def get_barcodes_file_list(filestem:str,num_files:int,padding:int)->list[str]:
    filenames = []
    for i in range(1,num_files+1):
        numstr = str(i)
        filenames.append(filestem+"0"*(padding-len(numstr))+numstr+".tsv")
    return filenames

def construct_dataframe(barcode_files:list[str])->pl.DataFrame:
    partial_dfs = [pl.scan_csv(bcd_file, separator="\t") for bcd_file in barcode_files]
    full_df = pl.concat(partial_dfs).group_by("Barcode").agg(pl.col("count").sum()).sort("count",descending=True).collect()
    full_df = full_df.with_row_index("rank",offset=1)
    return full_df

def threshold_ordmag(sorted_counts:list[pl.Int64],n_expected:int)->int:
    threshold = int(round(sorted_counts[int(n_expected*0.01)] / 10,0))
    return threshold

def apply_thresholding(counts_df:pl.DataFrame,min_readcount:int,ordmag_expected:int)->pl.DataFrame:
    otsu_threshold = int(round(threshold_otsu(counts_df.get_column("count").to_numpy()),0))
    ordmag_threshold = threshold_ordmag(counts_df.get_column("count").to_numpy(),ordmag_expected)
    counts_df = counts_df.with_columns(
        passes_otsu=pl.col("count") >= otsu_threshold,
        passes_ordmag=pl.col("count") >= ordmag_threshold,
        passes_readcount=pl.col("count") >= min_readcount
    )
    counts_df = counts_df.with_columns(
        QC_pass = pl.col("passes_otsu") & pl.col("passes_ordmag") & pl.col("passes_readcount")
    )       
        
    return counts_df

def plot_filter_criterion(counts_df:pl.DataFrame,filter_col:str,title:str,ax)->None:
    passed_df = counts_df.filter(pl.col(filter_col))
    failed_df = counts_df.filter(~pl.col(filter_col))
    ax.scatter(x=passed_df.get_column("rank"),y=passed_df.get_column("count"),c="blue",label="Passed",s=1)
    ax.scatter(x=failed_df.get_column("rank"),y=failed_df.get_column("count"),c="gray",label="Failed",s=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode Rank")
    ax.set_ylabel("Readcount")
    title += f" (Found: {passed_df.shape[0]})"
    ax.set_title(title)

def plot_data(counts_df:pl.DataFrame,stats_dir:str,slx:str,sample:str,cellranger_barcodes:set[str],expected_cells:int)->None:
    statsdir = Path(stats_dir)
    if not statsdir.is_dir():
        statsdir.mkdir()
    #write full file to disk
    counts_df.write_csv(statsdir/f"{slx}.{sample}.barcode_stats.tsv",separator="\t")
    #plot barcode rank plots
    if cellranger_barcodes:
        fig,axs = plt.subplots(3,2,figsize = (16,24))
    else:
        fig,axs = plt.subplots(2,2,figsize=(16,16))
    plot_filter_criterion(counts_df,"passes_readcount","Minimum Readcount",axs[0,0])
    plot_filter_criterion(counts_df,"passes_ordmag",f"OrdMag Thresholding (Exp:{expected_cells})",axs[0,1])
    plot_filter_criterion(counts_df,"passes_otsu","Otsu Tresholding",axs[1,0])
    plot_filter_criterion(counts_df,"QC_pass","All Filters",axs[1,1])
    if cellranger_barcodes:
        plot_filter_criterion(counts_df,"in_cellranger","Cellranger calls",axs[2,0])
        venn2(
            [set(counts_df.filter(pl.col("QC_pass")).get_column("RNA")),cellranger_barcodes],
            ("DNA-based","Cellranger"),
            ax=axs[2,1]
        )
    plt.tight_layout()
    plt.savefig(statsdir/f"{slx}.{sample}.barcode_ranks.png",dpi=600)




def produce_cell_call_file(counts_df:pl.DataFrame,outdir:str,slx:str,sample:str,flowcell:str,lane:str)->None:
    output_dir = Path(outdir)
    if not output_dir.is_dir():
        output_dir.mkdir()

    file = output_dir / f"{slx}.{sample}.{flowcell}.s_{lane}.selected.txt"
    counts_df.filter("QC_pass").select("Barcode").write_csv(file,include_header=False,separator="\t")

def main()->None:
    cellranger_barcodes = set()
    parser = get_parser()
    args = parser.parse_args()
    samples = list(zip(args.flowcell.split(" "),args.lane.split(" ")))
    print(samples)
    barcode_files = []
    for sample in samples:
        filestem = args.barcode_filestem + f"{sample[0]}.s_{sample[1]}."
        print(filestem)
        barcode_files += get_barcodes_file_list(filestem,args.num_files,args.padding)

    barcode_counts = construct_dataframe(barcode_files)

    barcode_counts = apply_thresholding(
        counts_df=barcode_counts,
        min_readcount=args.min_readcount,
        ordmag_expected=args.ordmag_expected
        )  
    
    plot_data(
        counts_df=barcode_counts,
        stats_dir=args.stats_dir,
        slx=args.slx,
        sample=args.sample,
        expected_cells=args.ordmag_expected,
        cellranger_barcodes=cellranger_barcodes
    )

    for sample in samples:
        produce_cell_call_file(
            counts_df=barcode_counts,
            outdir=args.output_dir,
            slx=args.slx,
            sample=args.sample,
            flowcell=sample[0],
            lane=sample[1]
        )



if __name__ == "__main__":
    main()
