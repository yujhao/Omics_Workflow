#!/usr/bin/env python3
# coding: utf-8

import os
from tools.load import load_diffgene_data, load_stringdb
from tools.utils import get_network, getPPIdatabase

outdir = "ppi_result"
os.makedirs(outdir, exist_ok=True)

diffgene = load_diffgene_data(
    "diffgene.xls",
    noft=None,
    gene_col="gene",
    fold_change_col="FoldChange",
    regulation_col="Regulation",
    up_regulation_value="Up",
    down_regulation_value="Down"
)

stringdb_filepath = getPPIdatabase(
    9606,
    "/home/liuchenglong/BuildPPIdatabase/ppi_databaseIdx.csv",
    "DatabasePath",
    "stringdbID"
)

stringdb = load_stringdb(stringdb_filepath)

network_df_top = get_network(
    diffgene,
    stringdb,
    300,
    "gene",
    "Regulation"
)

network_df = get_network(
    diffgene,
    stringdb,
    None,
    "gene",
    "Regulation"
)

network_df.to_csv(
    os.path.join(outdir, "result.ppi_network.tsv"),
    sep="\t",
    index=False
)

network_df_top.to_csv(
    os.path.join(outdir, "result.ppi_network_top300.tsv"),
    sep="\t",
    index=False
)

os.system(
    "bash script/network_3d.sh ppi_result result.ppi_network_top300.tsv"
)