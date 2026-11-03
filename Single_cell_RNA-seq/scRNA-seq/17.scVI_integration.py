#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# scVI (single-cell Variational Inference) integration and clustering
# Learns a batch-corrected latent representation across samples.
# Requirement: pip install scvi-tools scanpy
# ==============================================================================

import os

import numpy as np
import scanpy as sc
import scvi
from scvi.model.utils import mde

# ========================== Configuration parameters ==========================

OUTPUT_DIR = "1.scVI_integration"
INPUT_FILE = "data/data_raw.h5ad"       # raw count matrix in h5ad format
OUTPUT_FILE = "data/data_scvi.h5ad"

BATCH_KEY = "sampleid"                  # metadata column holding the sample id
N_TOP_GENES = 3000                      # number of highly variable genes
N_LATENT = 20                           # dimension of the latent space
RESOLUTION = 0.8                        # Leiden resolution
SEED = 0

# ========================== Setup ==========================

os.makedirs(OUTPUT_DIR, exist_ok=True)

scvi.settings.seed = SEED
scvi.settings.num_threads = 12
scvi.settings.device = "cuda"           # use "cpu" if no GPU is available

# ========================== Data loading ==========================

adata = sc.read_h5ad(INPUT_FILE)
adata = adata.raw.to_adata()

# scVI requires raw counts, keep them in a dedicated layer
adata.layers["counts"] = adata.X.copy()

# ========================== Normalization and HVG selection ==========================

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, batch_key=BATCH_KEY)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# ========================== Train the scVI model ==========================

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH_KEY)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=N_LATENT)
vae.train()

# ========================== Latent space and clustering ==========================

adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=50, key_added="scvi_neighbors")
sc.tl.leiden(adata, resolution=RESOLUTION, random_state=SEED, n_iterations=-1,
             directed=False, key_added="clusters_scvi", neighbors_key="scvi_neighbors")

sc.pl.embedding(adata, basis="X_mde", color=[BATCH_KEY, "clusters_scvi"],
                frameon=False, ncols=1, save="_scvi.png")

# ========================== Save results ==========================

adata.write(OUTPUT_FILE)
print("scVI integration finished, results written to", OUTPUT_FILE)
