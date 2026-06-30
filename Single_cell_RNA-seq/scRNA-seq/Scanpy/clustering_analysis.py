#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第二步：降维聚类分析
输入：QC 后的数据
功能：PCA / PCA+Harmony + 聚类 + UMAP/tSNE
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import argparse
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
LOGGER = logging.getLogger(__name__)


def run_pca(
    adata,
    n_comps=50,
    n_pcs=30,
    max_value=10
):
    """
    运行 PCA 降维

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象（已 QC）
    n_comps : int
        PCA 组件数
    n_pcs : int
        用于后续分析的 PC 数
    max_value : float
        scaling 的最大值
    """
    LOGGER.info("开始 PCA 降维...")

    # 如果还没有 HVG，先选择
    if 'highly_variable' not in adata.var.columns:
        LOGGER.info("选择高变基因...")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

    # 只使用高变基因
    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # Scaling
    sc.pp.scale(adata_hvg, max_value=max_value)
    LOGGER.info(f"Scaling 完成，max_value={max_value}")

    # PCA
    sc.tl.pca(adata_hvg, n_comps=n_comps)
    LOGGER.info(f"PCA 完成，{n_comps} 个组件")

    # 将 PCA 结果复制回原 adata
    # 注意：varm 和 uns 不需要复制，因为只需要 obsm 中的 PCA 坐标用于后续分析
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

    return adata, n_pcs


def run_harmony(
    adata,
    batch_key="batchid",
    n_pcs=30
):
    """
    运行 Harmony 批次校正

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象（已完成 PCA）
    batch_key : str
        批次键名
    n_pcs : int
        用于后续分析的 PC 数
    """
    LOGGER.info(f"开始 Harmony 批次校正，batch_key={batch_key}...")

    try:
        import harmony_pytorch as harmony

        # 运行 Harmony
        sc.external.pp.harmony_integrate(
            adata,
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony'
        )

        LOGGER.info("Harmony 批次校正完成")
        use_rep = 'X_pca_harmony'

    except ImportError:
        LOGGER.warning("未安装 harmony_pytorch，使用普通 PCA")
        use_rep = 'X_pca'

    return adata, use_rep


def run_clustering(
    adata,
    n_pcs=30,
    n_neighbors=15,
    resolution=0.8,
    method="leiden",
    use_rep="X_pca"
):
    """
    运行聚类

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    n_pcs : int
        用于聚类的 PC 数
    n_neighbors : int
        邻居数
    resolution : float
        聚类分辨率
    method : str
        聚类方法 (leiden 或 louvain)
    use_rep : str
        用于构建邻域图的表示
    """
    LOGGER.info(f"开始聚类，方法={method}, resolution={resolution}...")

    # 构建邻域图
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep
    )
    LOGGER.info(f"邻域图构建完成，n_neighbors={n_neighbors}, n_pcs={n_pcs}")

    # 聚类
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
        # 如果 clusters 已存在，先删除
        if 'clusters' in adata.obs.columns:
            adata.obs.drop(columns=['clusters'], inplace=True)
        adata.obs.rename(columns={"leiden": "clusters"}, inplace=True)
    elif method == "louvain":
        sc.tl.louvain(adata, resolution=resolution)
        # 如果 clusters 已存在，先删除
        if 'clusters' in adata.obs.columns:
            adata.obs.drop(columns=['clusters'], inplace=True)
        adata.obs.rename(columns={"louvain": "clusters"}, inplace=True)

    # 确保 clusters 是 categorical 类型
    adata.obs['clusters'] = pd.Categorical(adata.obs['clusters'])

    n_clusters = adata.obs['clusters'].nunique()
    LOGGER.info(f"聚类完成，共 {n_clusters} 个 clusters")

    return adata


def run_embedding(
    adata,
    n_pcs=30,
    min_dist=0.5,
    use_rep="X_pca",
    embedding_method="umap"
):
    """
    运行降维可视化

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    n_pcs : int
        用于 tSNE 的 PC 数
    min_dist : float
        UMAP 的 min_dist 参数
    use_rep : str
        用于 tSNE 的表示
    embedding_method : str
        降维方法 (umap, tsne, both)
    """
    LOGGER.info(f"开始降维可视化，方法={embedding_method}...")

    if embedding_method in ["umap", "both"]:
        sc.tl.umap(adata, min_dist=min_dist)
        LOGGER.info(f"UMAP 完成，min_dist={min_dist}")

    if embedding_method in ["tsne", "both"]:
        sc.tl.tsne(adata, n_pcs=n_pcs, use_rep=use_rep)
        LOGGER.info("tSNE 完成")

    return adata


def find_marker_genes(
    adata,
    method="wilcoxon",
    n_genes=100
):
    """
    查找每个 cluster 的 marker 基因

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    method : str
        检测方法 (wilcoxon, t-test, logreg)
    n_genes : int
        每个 cluster 返回的基因数

    Returns:
    --------
    DataFrame
        marker 基因结果
    """
    LOGGER.info(f"查找 marker 基因，方法={method}...")

    # 运行差异分析
    sc.tl.rank_genes_groups(
        adata,
        groupby='clusters',
        method=method,
        n_genes=n_genes,
        use_raw=False
    )

    # 提取结果
    result = adata.uns['rank_genes_groups']

    # 转换为 DataFrame
    all_results = []
    for cluster in result['names'].dtype.names:
        cluster_data = pd.DataFrame({
            'cluster': cluster,
            'gene': result['names'][cluster],
            'scores': result['scores'][cluster],
            'logfoldchanges': result['logfoldchanges'][cluster],
            'pvals': result['pvals'][cluster],
            'pvals_adj': result['pvals_adj'][cluster],
        })
        all_results.append(cluster_data)

    marker_df = pd.concat(all_results, ignore_index=True)

    LOGGER.info(f"找到 {len(marker_df)} 个 marker 基因")

    return marker_df


def add_gene_annotation(df, anno_df, gene_col="gene"):
    """
    为基因添加注释

    Parameters:
    -----------
    df : DataFrame
        基因数据表
    anno_df : DataFrame
        基因注释表
    gene_col : str
        基因名列名

    Returns:
    --------
    DataFrame
        添加注释后的数据表
    """
    # 合并注释
    df_merged = df.merge(
        anno_df,
        left_on=gene_col,
        right_on="id",
        how="left"
    )

    # 删除重复的 id 列
    if "id" in df_merged.columns:
        df_merged = df_merged.drop(columns=["id"])

    # 填充缺失的注释为 '--'
    df_merged = df_merged.fillna('--')

    # 重新排列列顺序
    base_cols = ["cluster", "gene", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    anno_cols = ["ensembl_id", "gene_type", "gene_description", "TFs_Family", "GO_id", "GO_term", "KEGG_id", "KEGG_description"]

    # 保留存在的列
    ordered_cols = [col for col in base_cols + anno_cols if col in df_merged.columns]
    other_cols = [col for col in df_merged.columns if col not in ordered_cols]
    ordered_cols.extend(other_cols)

    df_merged = df_merged[ordered_cols]

    return df_merged


def plot_embeddings(
    adata,
    output_dir,
    color_by=["clusters", "sampleid"],
    embedding_method="umap"
):
    """
    绘制降维可视化图

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    output_dir : str
        输出目录
    color_by : list
        用于着色的列名
    embedding_method : str
        降维方法 (umap, tsne, both)
    """
    LOGGER.info("绘制降维可视化图...")

    os.makedirs(output_dir, exist_ok=True)

    # 设置 scanpy 图片保存目录
    sc.settings.figdir = output_dir

    for color in color_by:
        if color not in adata.obs.columns:
            LOGGER.warning(f"列 {color} 不存在，跳过")
            continue

        try:
            if embedding_method in ["umap", "both"]:
                # UMAP 图
                sc.pl.umap(
                    adata,
                    color=color,
                    title=f'UMAP - {color}',
                    save=f'_umap_{color}.png',
                    show=False
                )
                LOGGER.info(f"已保存 UMAP 图：umap_{color}.png")

            if embedding_method in ["tsne", "both"]:
                # tSNE 图
                sc.pl.tsne(
                    adata,
                    color=color,
                    title=f'tSNE - {color}',
                    save=f'_tsne_{color}.png',
                    show=False
                )
                LOGGER.info(f"已保存 tSNE 图：tsne_{color}.png")
        except Exception as e:
            LOGGER.warning(f"绘制 {color} 时出错：{e}，跳过")


def save_results(adata, marker_df, outdir, refgenome=None):
    """
    保存结果

    Parameters:
    -----------
    adata : AnnData
        分析完成的 AnnData 对象
    marker_df : DataFrame
        marker 基因结果
    outdir : str
        输出目录
    refgenome : str
        参考基因组目录路径（用于基因注释）
    """
    LOGGER.info(f"保存结果到 {outdir}...")

    os.makedirs(outdir, exist_ok=True)

    # 保存 h5ad
    adata.write_h5ad(os.path.join(outdir, "clustering_result.h5ad"))
    LOGGER.info("已保存 clustering_result.h5ad")

    # 保存 UMAP 坐标（格式：Barcode, UMAP_1, UMAP_2）
    if 'X_umap' in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm['X_umap'],
            columns=['UMAP_1', 'UMAP_2'],
            index=adata.obs_names
        )
        umap_df.index.name = 'Barcode'
        umap_df.to_csv(os.path.join(outdir, "umap_coordinates.csv"))
        LOGGER.info("已保存 umap_coordinates.csv")

    # 保存聚类结果（格式：Barcode, sampleid, group, batchid, clusters）
    cluster_result_df = pd.DataFrame(index=adata.obs_names)
    cluster_result_df.index.name = 'Barcode'

    # 添加可用的列
    for col in ['sampleid', 'batchid', 'clusters']:
        if col in adata.obs.columns:
            cluster_result_df[col] = adata.obs[col]

    # 如果有 group 列则添加，否则用 batchid 作为 group
    if 'group' in adata.obs.columns:
        cluster_result_df['group'] = adata.obs['group']
    elif 'batchid' in adata.obs.columns:
        cluster_result_df['group'] = adata.obs['batchid']

    # 重新排列列顺序
    desired_order = ['sampleid', 'group', 'batchid', 'clusters']
    cols = [col for col in desired_order if col in cluster_result_df.columns]
    cluster_result_df = cluster_result_df[cols]

    cluster_result_df.to_csv(os.path.join(outdir, "clusters_result.csv"))
    LOGGER.info("已保存 clusters_result.csv")

    # 添加基因注释并保存 marker 基因
    if refgenome and os.path.exists(os.path.join(refgenome, "gene_annotation.xls")):
        LOGGER.info(f"加载基因注释: {refgenome}/gene_annotation.xls")
        anno_df = pd.read_csv(os.path.join(refgenome, "gene_annotation.xls"), sep="\t")
        anno_df["id"] = anno_df["id"].astype(str)

        # 添加注释
        marker_df_annotated = add_gene_annotation(marker_df, anno_df)

        # 保存所有 marker 基因
        marker_df_annotated.to_csv(os.path.join(outdir, "marker_genes.csv"), index=False)
        LOGGER.info("已保存 marker_genes.csv (带基因注释)")

        # 保存每个 cluster 的 top 10 marker 基因
        top_markers = marker_df_annotated.groupby("cluster").head(10)
        top_markers.to_csv(os.path.join(outdir, "top_marker_genes.csv"), index=False)
        LOGGER.info("已保存 top_marker_genes.csv (带基因注释)")
    else:
        # 没有注释文件，保存原始数据
        marker_df.to_csv(os.path.join(outdir, "marker_genes.csv"), index=False)
        LOGGER.info("已保存 marker_genes.csv")

        # 保存每个 cluster 的 top 10 marker 基因
        top_markers = marker_df.groupby("cluster").head(10)
        top_markers.to_csv(os.path.join(outdir, "top_marker_genes.csv"), index=False)
        LOGGER.info("已保存 top_marker_genes.csv")

    LOGGER.info("所有结果保存完成")


def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(
        description="第二步：降维聚类分析（PCA / PCA+Harmony + 聚类 + UMAP/tSNE + Marker 基因）"
    )

    # 输入输出
    parser.add_argument("-i", "--input", required=True, help="输入文件路径 (QC 后的 .h5ad)")
    parser.add_argument("-o", "--outdir", required=True, help="输出目录")

    # PCA 参数
    parser.add_argument("--n_comps", type=int, default=50, help="PCA 组件数 (default: 50)")
    parser.add_argument("--n_pcs", type=int, default=30, help="用于后续分析的 PC 数 (default: 30)")
    parser.add_argument("--max_value", type=float, default=10.0, help="scaling 最大值 (default: 10)")

    # Harmony 参数
    parser.add_argument("--use_harmony", action="store_true", help="是否使用 Harmony 批次校正")
    parser.add_argument("--batch_key", type=str, default="batchid", help="批次键名 (default: batchid)")

    # 聚类参数
    parser.add_argument("--cluster_method", type=str, default="leiden", choices=["leiden", "louvain"],
                        help="聚类方法 (default: leiden)")
    parser.add_argument("--resolution", type=float, default=0.8, help="聚类分辨率 (default: 0.8)")
    parser.add_argument("--n_neighbors", type=int, default=15, help="邻居数 (default: 15)")

    # 可视化参数
    parser.add_argument("--embedding", type=str, default="umap", choices=["umap", "tsne", "both"],
                        help="降维可视化方法 (default: umap)")
    parser.add_argument("--min_dist", type=float, default=0.5, help="UMAP min_dist (default: 0.5)")

    # Marker 基因参数
    parser.add_argument("--marker_method", type=str, default="wilcoxon",
                        choices=["wilcoxon", "t-test", "logreg"],
                        help="Marker 基因检测方法 (default: wilcoxon)")
    parser.add_argument("--n_genes", type=int, default=100, help="每个 cluster 返回的基因数 (default: 100)")
    parser.add_argument("--refgenome", type=str, default=None, help="参考基因组目录路径（用于基因注释）")

    args = parser.parse_args()

    # 读取数据
    LOGGER.info(f"读取数据：{args.input}")
    adata = sc.read_h5ad(args.input)

    # 1. PCA 降维
    adata, n_pcs = run_pca(
        adata,
        n_comps=args.n_comps,
        n_pcs=args.n_pcs,
        max_value=args.max_value
    )

    # 2. Harmony 批次校正（可选）
    if args.use_harmony:
        adata, use_rep = run_harmony(
            adata,
            batch_key=args.batch_key,
            n_pcs=args.n_pcs
        )
    else:
        use_rep = "X_pca"
        LOGGER.info("不使用 Harmony，使用普通 PCA")

    # 3. 聚类
    adata = run_clustering(
        adata,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        method=args.cluster_method,
        use_rep=use_rep
    )

    # 4. 降维可视化
    adata = run_embedding(
        adata,
        n_pcs=args.n_pcs,
        min_dist=args.min_dist,
        use_rep=use_rep,
        embedding_method=args.embedding
    )

    # 5. 绘制可视化图
    color_by = ["clusters"]
    if "sampleid" in adata.obs.columns:
        color_by.append("sampleid")
    if "batchid" in adata.obs.columns and "batchid" != "sampleid":
        color_by.append("batchid")

    plot_embeddings(
        adata,
        os.path.join(args.outdir, "embeddings"),
        color_by=color_by,
        embedding_method=args.embedding
    )

    # 6. 查找 marker 基因
    marker_df = find_marker_genes(
        adata,
        method=args.marker_method,
        n_genes=args.n_genes
    )

    # 7. 保存结果
    save_results(adata, marker_df, args.outdir, refgenome=args.refgenome)

    LOGGER.info("分析完成！")


if __name__ == "__main__":
    main()
