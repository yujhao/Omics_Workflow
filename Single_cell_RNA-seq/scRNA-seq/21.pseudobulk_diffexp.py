#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Pseudobulk differential expression analysis
# Cells of the same sample and cell type are summed into pseudobulk profiles,
# which are then tested for differential expression between two conditions with
# DESeq2. Aggregating cells restores biological replicates and removes the
# pseudo-replication problem of cell level tests.
# Requirement: pip install decoupler pydeseq2 scanpy anndata
# ==============================================================================

import os

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DefaultInference, DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.sparse import issparse

# ========================== Configuration parameters ==========================

OUTPUT_DIR = "21.pseudobulk_diffexp"
INPUT_FILE = "data/cluster.h5ad"        # h5ad after clustering

SAMPLE_COL = "sampleid"                 # biological replicate column
GROUP_COL = "clusters"                  # cell type column, None to aggregate by sample only
CONDITION_COL = "group"                 # sample level column holding the two conditions
TREATMENT = "Treatment"                 # level of CONDITION_COL used as treatment
CONTROL = "Control"                     # level of CONDITION_COL used as control

MIN_CELLS = 10                          # min cells per sample x group combination
MIN_COUNTS = 1000                       # min total counts per pseudobulk sample
TARGET_SUM = 1e4                        # target sum used to recompute fold changes

FIT_TYPE = "mean"                       # parametric, local or mean
N_CPUS = 8
ALPHA = 0.05

# ========================== Helper functions ==========================


def standardize_expression_layers(adata):
    """Make sure adata.X holds raw counts, which DESeq2 requires."""
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    values = adata.X.data if issparse(adata.X) else np.asarray(adata.X).ravel()
    if not np.all(np.mod(values, 1) == 0):
        raise ValueError("Non-integer values detected: adata.X is not a raw count matrix.")

    return adata


def normalize_total_log1p(pdata, target_sum):
    """Normalize each pseudobulk sample to target_sum, then apply log1p."""
    X = pdata.X.toarray() if issparse(pdata.X) else np.array(pdata.X, dtype=float)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    return np.log1p(X / totals * target_sum)


def compute_seurat_logfc(X_group, X_rest):
    """Seurat style log2 fold change: ratio of mean exponentiated expression."""
    mean_group = np.mean(np.expm1(X_group), axis=0) + 1.0 / X_group.shape[0]
    mean_rest = np.mean(np.expm1(X_rest), axis=0) + 1.0 / X_rest.shape[0]

    fold_changes = mean_group / mean_rest
    log_fc = np.log2(fold_changes)
    log_fc[np.isnan(log_fc)] = 0

    return log_fc, fold_changes


# ========================== Initialization ==========================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========================== Data loading ==========================

adata = standardize_expression_layers(sc.read_h5ad(INPUT_FILE))

groups_col = None if GROUP_COL in (None, "None") else GROUP_COL
for col in [SAMPLE_COL] + ([groups_col] if groups_col else []):
    if col not in adata.obs:
        raise ValueError(f"Column '{col}' not found in adata.obs")

# ========================== Pseudobulk aggregation ==========================

# Every sample x group combination becomes one pseudobulk sample
pdata = dc.pp.pseudobulk(adata, sample_col=SAMPLE_COL, groups_col=groups_col, mode="sum")

# Remove pseudobulk samples below the quality thresholds
dc.pp.filter_samples(pdata, min_cells=MIN_CELLS, min_counts=MIN_COUNTS)
print(f"Pseudobulk samples retained: {pdata.n_obs}")

# ========================== DESeq2 ==========================

# DESeq2 takes raw counts directly, no normalization here
inference = DefaultInference(n_cpus=N_CPUS)

dds = DeseqDataSet(
    adata=pdata,
    design=f"~{CONDITION_COL}",
    refit_cooks=True,
    inference=inference,
    size_factors_fit_type="ratio",   # equivalent to estimateSizeFactors in R
    fit_type=FIT_TYPE,
    min_replicates=7,                # equivalent to minReplicatesForReplace
    beta_tol=1e-8,
    min_mu=0.5,
    min_disp=1e-8,
    max_disp=10.0,
)
dds.deseq2(fit_type=FIT_TYPE)

stat_res = DeseqStats(
    dds=dds,
    contrast=[CONDITION_COL, TREATMENT, CONTROL],
    inference=inference,
    alpha=ALPHA,
    cooks_filter=True,
    independent_filter=True,
)
stat_res.summary()

results_df = stat_res.results_df
results_df.index.name = "gene"
results_df = results_df.reset_index(drop=False)

# Use the same column names as the regular wilcoxon results
results_df = results_df.rename(columns={"pvalue": "p-value", "padj": "q-value"})

# ========================== Recompute fold changes and detection rates ==========================

# Normalize the pseudobulk matrix, DESeq2 size factors are not used for logFC
norm_log = normalize_total_log1p(pdata, TARGET_SUM)

group_mask = (pdata.obs[CONDITION_COL] == TREATMENT).values
rest_mask = (pdata.obs[CONDITION_COL] == CONTROL).values

log_fc, fold_changes = compute_seurat_logfc(norm_log[group_mask], norm_log[rest_mask])

# Align the statistics to the gene order of the results table
gene_indices = [pdata.var_names.get_loc(gene) for gene in results_df["gene"]]

results_df["FoldChange"] = fold_changes[gene_indices]
results_df["log2FoldChange"] = log_fc[gene_indices]

# pct.1 / pct.2: fraction of pseudobulk samples with non-zero expression
X_counts = pdata.X.toarray() if issparse(pdata.X) else np.asarray(pdata.X)
results_df["pct.1"] = np.round((X_counts[group_mask] > 0).mean(axis=0)[gene_indices], 3)
results_df["pct.2"] = np.round((X_counts[rest_mask] > 0).mean(axis=0)[gene_indices], 3)

results_df = results_df.drop(columns=["lfcSE", "stat"])

# Statistics that could not be computed are reported as 1
results_df["p-value"] = results_df["p-value"].fillna(1.0)
results_df["q-value"] = results_df["q-value"].fillna(1.0)

column_order = ["gene", "p-value", "log2FoldChange", "pct.1", "pct.2",
                "q-value", "FoldChange", "baseMean"]
results_df = results_df[column_order]

results_df.to_csv(os.path.join(OUTPUT_DIR, "pseudobulk_all_diffexp_genes.xls"),
                  sep="\t", index=False)

print(f"Pseudobulk differential expression finished, results written to {OUTPUT_DIR}")