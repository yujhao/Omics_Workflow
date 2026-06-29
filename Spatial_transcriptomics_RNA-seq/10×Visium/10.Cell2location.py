#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import cell2location

rcParams["pdf.fonttype"] = 42

# =========================
# 1. 路径和常用参数
# =========================

sc_file = "./sc.h5ad"
st_file = "./adata_sp.h5ad"
outdir = "./result"

labels_key = "new_celltype"
batch_key = "sampleid"

n_cells_per_location = 10
detection_alpha = 20

os.makedirs(outdir, exist_ok=True)

# =========================
# 2. 读数据
# =========================

adata_ref = sc.read_h5ad(sc_file)
adata_vis = sc.read_h5ad(st_file)

# =========================
# 3. 去掉线粒体基因
# =========================

adata_ref.obsm["mt"] = adata_ref[:, adata_ref.var["mt"].values].X.toarray()
adata_ref = adata_ref[:, ~adata_ref.var["mt"].values].copy()

adata_vis.obsm["mt"] = adata_vis[:, adata_vis.var["mt"].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var["mt"].values].copy()

# =========================
# 4. 单细胞参考过滤基因
# =========================

from cell2location.utils.filtering import filter_genes

selected = filter_genes(
    adata_ref,
    cell_count_cutoff=5,
    cell_percentage_cutoff2=0.03,
    nonz_mean_cutoff=1.12
)

adata_ref = adata_ref[:, selected].copy()

# =========================
# 5. 训练单细胞 reference model
# =========================

from cell2location.models import RegressionModel

RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key=batch_key,
    labels_key=labels_key,
    categorical_covariate_keys=None
)

mod = RegressionModel(adata_ref)
mod.view_anndata_setup()

mod.train(
    max_epochs=150,
    batch_size=1000,
    accelerator="cpu"
)

mod.plot_history(20)
plt.savefig(os.path.join(outdir, "01.reference_history.pdf"))
plt.close()

mod.plot_QC()
plt.savefig(os.path.join(outdir, "02.reference_QC.pdf"))
plt.close()

adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": "cpu"}
)

mod.save(os.path.join(outdir, "reference_model"), overwrite=True)
adata_ref.write(os.path.join(outdir, "sc_reference.h5ad"))

# =========================
# 6. 提取 cell type signatures
# =========================

if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()

inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
inf_aver.to_csv(os.path.join(outdir, "cell_type_signatures.csv"))

# =========================
# 7. 共同基因
# =========================

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# =========================
# 8. 空间映射
# =========================

cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis,
    batch_key=batch_key
)

mod = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    N_cells_per_location=n_cells_per_location,
    detection_alpha=detection_alpha
)

mod.view_anndata_setup()

mod.train(
    max_epochs=30000,
    batch_size=None,
    train_size=1,
    accelerator="cpu"
)

mod.plot_history(1000)
plt.legend(labels=["full data training"])
plt.savefig(os.path.join(outdir, "03.mapping_history.pdf"))
plt.close()

adata_vis = mod.export_posterior(
    adata_vis,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": mod.adata.n_obs,
        "accelerator": "cpu"
    }
)

mod.save(os.path.join(outdir, "mapping_model"), overwrite=True)
adata_vis.write(os.path.join(outdir, "st_cell2location.h5ad"))

mod.plot_QC()
plt.savefig(os.path.join(outdir, "04.mapping_QC.png"))
plt.close()

# =========================
# 9. 导出 abundance
# =========================

abundance = adata_vis.obsm["q05_cell_abundance_w_sf"].copy()
abundance.columns = adata_vis.uns["mod"]["factor_names"]
abundance.to_csv(os.path.join(outdir, f"{sampleid}.q05_cell_abundance.csv"))

adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = abundance

# =========================
# 10. 按样本取一个 slide 出图
# =========================

from cell2location.utils import select_slide

slide = select_slide(adata_vis, sampleid, batch_key=batch_key)

with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [4.5, 5]}):
    sc.pl.spatial(
        slide,
        cmap="magma",
        color=list(adata_vis.uns["mod"]["factor_names"])[:9],
        ncols=4,
        size=1.3,
        img_key="hires",
        vmin=0,
        vmax="p99.2",
        show=False
    )

plt.savefig(os.path.join(outdir, "05.spatial_celltype.pdf"))
plt.close()

# =========================
# 11. 多细胞类型合图
# =========================

from cell2location.plt import plot_spatial

fig = plot_spatial(
    adata=slide,
    color=list(adata_vis.uns["mod"]["factor_names"])[:9],
    labels=list(adata_vis.uns["mod"]["factor_names"])[:9],
    show_img=True,
    img_key="hires",
    style="fast",
    max_color_quantile=0.992,
    circle_diameter=6,
    colorbar_position="right"
)

plt.savefig(os.path.join(outdir, "06.spatial_all_celltypes.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(outdir, "06.spatial_all_celltypes.png"), dpi=300, bbox_inches="tight")
plt.close()