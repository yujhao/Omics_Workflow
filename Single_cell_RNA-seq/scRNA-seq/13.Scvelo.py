#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
scVelo RNA 速度分析脚本
用于单细胞 RNA-seq 数据的速度分析和可视化

前置步骤：
    velocyto run10x Cellranger/sample/ genes.gtf
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvelo as scv

# ========================== Configuration parameters ==========================

# Input file paths
LOOM_FILE = "utils/data/test.loom"
CELLS_FILE = "utils/data/subset_cells.txt"
UMAP_FILE = "utils/data/umap_coords.csv"
CELLTYPES_FILE = "utils/data/cell_types.csv"

# Output directory
OUTPUT_DIR = Path("13.Scvelo")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Analysis parameters
MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
N_JOBS = 10

# ========================== Data loading ==========================

# 1. Read loom file
vlm = sc.read_loom(LOOM_FILE, X_name='spliced')

# 2. Cell subset filtering
subset_cells = pd.read_csv(CELLS_FILE, header=None)[0].tolist()
vlm_subset = vlm[subset_cells, :].copy()

# 3. Import Seurat UMAP
umap_coords = pd.read_csv(UMAP_FILE, index_col=0)
vlm_subset.obsm['X_umap'] = umap_coords.loc[vlm_subset.obs_names].values

# 4. Import cell types
cell_types = pd.read_csv(CELLTYPES_FILE, index_col=0)
vlm_subset.obs['CellType'] = cell_types.loc[vlm_subset.obs_names, 'CellType']
vlm_subset.obs["CellType"] = vlm_subset.obs["CellType"].astype("category")

# ========================== Data preprocessing ==========================

scv.pp.filter_and_normalize(vlm_subset, min_shared_counts=MIN_SHARED_COUNTS, n_top_genes=N_TOP_GENES)
scv.pp.moments(vlm_subset)

# Plot spliced/unspliced proportion plot
scv.pl.proportions(
    vlm_subset,
    groupby="CellType",
    save=OUTPUT_DIR / "proportions_spliced_unspliced_counts_groupby.pdf",
    figsize=(12, 5),
    fontsize=10
)

# ========================== RNA Velocity calculation ==========================

# Compute dynamic velocity
scv.tl.recover_dynamics(vlm_subset, n_jobs=N_JOBS)
# Dynamic model
scv.tl.velocity(vlm_subset, vkey='dynamical_velocity', mode='dynamical')
scv.tl.velocity_graph(vlm_subset, vkey='dynamical_velocity', n_jobs=N_JOBS)

# ========================== RNA Velocity visualization ==========================

# 1. Stream plot
scv.pl.velocity_embedding_stream(
    vlm_subset,
    vkey='dynamical_velocity',
    basis='umap',
    color='CellType',
    legend_loc='right',
    figsize=(8, 6),
    dpi=150,
    show=False
)
plt.savefig(OUTPUT_DIR / "dynamical_velocity_groupby_CellType.svg", bbox_inches="tight", dpi=300)

# 2. Arrow plot
scv.pl.velocity_embedding(
    vlm_subset,
    basis='umap',
    vkey='dynamical_velocity',
    alpha=1,
    arrow_length=6,
    arrow_size=6,
    color='CellType',
    legend_fontsize=10,
    legend_loc='right margin',
    save=None
)
plt.savefig(OUTPUT_DIR / "cell_level_dynamical_velocity_groupby_CellType.svg", bbox_inches="tight", dpi=300)

# 3. Grid plot
scv.pl.velocity_embedding_grid(
    vlm_subset,
    basis='umap',
    vkey='dynamical_velocity',
    color='CellType',
    alpha=1,
    size=30,
    scale=0.25,
    title=None,
    legend_loc='right margin',
    legend_fontsize=10,
    save=None
)
plt.savefig(OUTPUT_DIR / "scvelo_embedding_grid_groupby_CellType.svg", bbox_inches="tight", dpi=300)

# ========================== Latent Time analysis ==========================

# Compute latent time
scv.tl.recover_latent_time(vlm_subset, vkey='dynamical_velocity')

# Plot latent time
scv.pl.scatter(
    vlm_subset,
    basis='umap',
    color='latent_time',
    fontsize=24,
    size=100,
    title=None,
    color_map='gnuplot',
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    save=None
)
plt.savefig(OUTPUT_DIR / "latent_time_by_dynamical.pdf", bbox_inches="tight", dpi=300)

# ========================== Gene analysis ==========================

# Get top genes
top_genes = vlm_subset.var['fit_likelihood'].sort_values(ascending=False).index[:20]

# Heatmap
scv.pl.heatmap(
    vlm_subset,
    var_names=top_genes,
    sortby='latent_time',
    col_color='CellType',
    n_convolve=300,
    save=None
)
plt.savefig(OUTPUT_DIR / "latent_time_heatmap.pdf", bbox_inches="tight", dpi=300)

# Top likelihood gene scatter plot
scv.pl.scatter(
    vlm_subset,
    basis=top_genes[:1],
    use_raw=False,
    color='CellType',
    frameon=False,
    show=False,
    legend_loc='right margin'
)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Top-likelihood_genes_groupby.pdf", dpi=300)

# Top likelihood gene time-trajectory plot
scv.pl.scatter(
    vlm_subset,
    x='latent_time',
    y=top_genes[:1],
    use_raw=False,
    color='CellType',
    frameon=False,
    legend_loc='right margin'
)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Top-likelihood_genes_alonetime_groupby.pdf", dpi=300)

# ========================== Gene ranking analysis ==========================

# Dynamic gene ranking
scv.tl.rank_dynamical_genes(vlm_subset, groupby='CellType')
df_dynamical = scv.get_df(vlm_subset, 'rank_dynamical_genes/names')
df_dynamical.to_excel(OUTPUT_DIR / "rank_dynamical_genes.xlsx", index=False)

# Velocity gene ranking
scv.tl.rank_velocity_genes(
    vlm_subset,
    vkey='dynamical_velocity',
    groupby='CellType',
    min_corr=0.3
)
df_velocity = pd.DataFrame(vlm_subset.uns['rank_velocity_genes']['names'])
df_velocity.to_excel(OUTPUT_DIR / "rank_velocity_genes.xlsx", index=False)

print(f"文件已成功保存到：{OUTPUT_DIR}")
print(f"1. 动态基因排名：{OUTPUT_DIR / 'rank_dynamical_genes.xlsx'}")
print(f"2. 速度基因排名：{OUTPUT_DIR / 'rank_velocity_genes.xlsx'}")

# ========================== Velocity confidence analysis ==========================

scv.tl.velocity_confidence(vlm_subset, vkey='dynamical_velocity')
keys = ['dynamical_velocity_length', 'dynamical_velocity_confidence']

scv.pl.scatter(
    vlm_subset,
    c=keys,
    cmap="coolwarm",
    basis='umap',
    perc=[5, 95],
    show=False
)
plt.savefig(OUTPUT_DIR / "Speed_and_coherence.pdf", dpi=300, bbox_inches="tight")
