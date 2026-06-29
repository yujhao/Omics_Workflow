#!/usr/bin/env python
# coding: utf-8
#from pathlib import Path
import argparse
#import re
# import cellphonedb as scv
#import pandas as pd
#import scanpy as sc
#import numpy as np
#import anndata2ri
#import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
# import difflib
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.pyplot as plt


cellphonedb = argparse.ArgumentParser(description='cellphonedb')
cellphonedb.add_argument('--counts', type=str, default = None)
cellphonedb.add_argument('--metadata', type=str, default = None)
cellphonedb.add_argument('--database', type=str, default = None)
cellphonedb.add_argument('--pvalue', type= float, default =0.05)
cellphonedb.add_argument('--threshold', type=float, default =0.1)
cellphonedb.add_argument('--threads', type=int, default =10 )
cellphonedb.add_argument('--iterations', type=int, default =1000 )
cellphonedb.add_argument('--microenvs', type=str,default=None)
cellphonedb.add_argument('--genetype', type=str, default = "gene_name")

args = cellphonedb.parse_args()
################################################解析参数##########################################################
cpdb_file_path = args.database
meta_file_path = args.metadata
counts_file_path = args.counts
pvalue=args.pvalue
out_path = "out"
threshold = args.threshold
genetype = args.genetype
iterations = args.iterations
if args.microenvs == "" :
  microenvs_file_path = None
else :
  microenvs_file_path = args.microenvs
###############################################主函数运行##########################################################
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
    counts_data = genetype,                     # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    iterations = iterations,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 4,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = ""                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )