# mres-data
The used code and data for the MRes dissertation "The Limits and Promise of Computational Integration: An Evaluation of CNV-Based Matching for Genomic and Transcriptomic Data at Single-Cell Resolution."
## scRNA-seq preprocessing
The snakemake workflow used for scRNA-seq pre-processing and CNV inference is in 'scRNA_cv/'
## Cell matching
All code for cell matching, both for cell to cell assessment and cell line to cell line assessment is in the 'cell_matching/' directory. The obtained accuracies can be found in 'postprocessing/accuracies/'
## Statistics and sumaries
The code used to summarise and plot the results and all statistical test on the accuracies is in 'postprocessing/'. All summaries and statistical test results are in 'postprocessing/statistics_and_summaries/', all figures are in 'postprocessing/figures/'