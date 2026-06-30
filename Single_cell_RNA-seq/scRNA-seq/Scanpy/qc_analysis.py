#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第一步：质量控制（QC）
输入：原始数据（10x 输出目录或多个样本）
功能：QC 过滤 + 双细胞去除
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


def load_gene_list(path):
    """
    从 GMT 文件加载基因列表
    """
    genes = []
    with open(path, "r") as f:
        for line in f:
            genes.extend(line.strip().split())

    genes = [g for g in genes if g and g != "NA"]
    return np.unique(genes)


def load_samples(root_dir, sample_dirs):
    """
    读取多个样本的 10x 数据并合并

    Parameters:
    -----------
    root_dir : str
        样本数据的根目录
    sample_dirs : list
        样本目录列表

    Returns:
    --------
    AnnData
        合并后的 AnnData 对象
    """
    LOGGER.info(f"从 {root_dir} 读取 {len(sample_dirs)} 个样本...")

    adata_list = []
    for sample in sample_dirs:
        sample_path = os.path.join(root_dir, sample)

        # 尝试读取 h5 文件
        h5_path = os.path.join(sample_path, "filtered_feature_bc_matrix.h5")
        mtx_path = os.path.join(sample_path, "filtered_feature_bc_matrix")

        if os.path.exists(h5_path):
            LOGGER.info(f"读取样本 {sample}: {h5_path}")
            adata = sc.read_10x_h5(h5_path)
        elif os.path.exists(mtx_path):
            LOGGER.info(f"读取样本 {sample}: {mtx_path}")
            adata = sc.read_10x_mtx(mtx_path, var_names='gene_symbols', cache=False)
        else:
            raise FileNotFoundError(f"找不到样本 {sample} 的数据文件")

        # 添加样本信息
        adata.obs['sampleid'] = sample
        adata.obs['batchid'] = sample  # 批次信息与样本相同

        # 确保 var_names 唯一
        adata.var_names_make_unique()

        adata_list.append(adata)

    # 合并所有样本
    LOGGER.info(f"合并 {len(adata_list)} 个样本...")
    adata_merged = ad.concat(adata_list, join='outer', merge='same')

    # 确保 obs_names 唯一
    adata_merged.obs_names_make_unique()

    LOGGER.info(f"合并完成: {adata_merged.n_obs} 个细胞, {adata_merged.n_vars} 个基因")

    return adata_merged


def calculate_qc_metrics(
    adata,
    refgenome,
    mt_genelist_file="MT_genelist.gmt",
    hb_genelist_file="HB_genelist.gmt"
):
    """
    计算质控指标 - 参考 qc.py 的实现

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    refgenome : str
        参考基因组目录路径
    mt_genelist_file : str
        线粒体基因列表文件名
    hb_genelist_file : str
        血红蛋白基因列表文件名
    """
    LOGGER.info("开始计算 QC 指标...")

    # 加载基因列表
    mt_genes = load_gene_list(os.path.join(refgenome, mt_genelist_file))
    hb_genes = load_gene_list(os.path.join(refgenome, hb_genelist_file))

    LOGGER.info(f"线粒体基因数: {len(mt_genes)}")
    LOGGER.info(f"血红蛋白基因数: {len(hb_genes)}")

    # 标记基因（使用 isin 而不是 startswith）
    adata.var["mt"] = adata.var_names.isin(mt_genes)
    adata.var["hb"] = adata.var_names.isin(hb_genes)

    # 计算 QC 指标
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "hb"],
        inplace=True
    )

    # 统一命名（Seurat 风格）
    adata.obs["nFeature_RNA"] = adata.obs["n_genes_by_counts"]
    adata.obs["nCount_RNA"] = adata.obs["total_counts"]
    adata.obs["percent.mito"] = adata.obs["pct_counts_mt"] / 100
    adata.obs["percent.HB"] = adata.obs["pct_counts_hb"] / 100

    LOGGER.info(f"QC 指标计算完成，共有 {adata.n_obs} 个细胞，{adata.n_vars} 个基因")

    return adata


def filter_cells(
    adata,
    mingene=200,
    maxgene=7500,
    minumi=1000,
    maxumi=50000,
    mtfilter=0.2,
    hbfilter=0.05
):
    """
    过滤细胞 - 参考 qc.py 的正确实现

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    mingene : int
        每个细胞最少基因数
    maxgene : int
        每个细胞最多基因数
    minumi : int
        每个细胞最少 UMI 数
    maxumi : int
        每个细胞最多 UMI 数
    mtfilter : float
        线粒体基因比例上限
    hbfilter : float
        血红蛋白基因比例上限
    """
    LOGGER.info("开始过滤细胞...")

    n_before = adata.n_obs

    # 过滤条件
    filter_mask = (
        (adata.obs["nFeature_RNA"] >= mingene) &
        (adata.obs["nFeature_RNA"] <= maxgene) &
        (adata.obs["nCount_RNA"] >= minumi) &
        (adata.obs["nCount_RNA"] <= maxumi) &
        (adata.obs["percent.mito"] <= mtfilter) &
        (adata.obs["percent.HB"] <= hbfilter)
    )

    # 应用过滤
    adata = adata[filter_mask].copy()

    n_after = adata.n_obs
    LOGGER.info(f"过滤完成：{n_before} -> {n_after} 个细胞 (移除 {n_before - n_after} 个)")

    return adata


def remove_doublets(adata, random_state=2025):
    """
    去除双细胞

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    random_state : int
        随机种子
    """
    LOGGER.info("开始双细胞检测...")

    try:
        import doubletdetection

        n_before = adata.n_obs

        # 保存原始计数
        adata.layers["counts"] = adata.X.copy()

        # 运行双细胞检测
        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            random_state=random_state
        )
        doublets = clf.fit(adata.X).predict()

        # 去除双细胞
        adata = adata[doublets == 0].copy()

        n_after = adata.n_obs
        LOGGER.info(f"双细胞去除完成：{n_before} -> {n_after} 个细胞 (移除 {n_before - n_after} 个)")

    except ImportError:
        LOGGER.warning("未安装 doubletdetection，跳过双细胞去除")

    return adata


def plot_qc_violin(adata, output_dir, title_suffix, group_col="sampleid"):
    """
    绘制 QC 小提琴图

    Parameters:
    -----------
    adata : AnnData
        输入的 AnnData 对象
    output_dir : str
        输出目录
    title_suffix : str
        标题后缀（如 "beforeQC" 或 "afterQC"）
    group_col : str
        分组列名
    """
    LOGGER.info(f"绘制 QC 小提琴图 ({title_suffix})...")

    os.makedirs(output_dir, exist_ok=True)

    metrics = ["nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"]

    # 检查分组列是否存在
    if group_col not in adata.obs.columns:
        LOGGER.warning(f"列 {group_col} 不存在，使用所有数据绘制")
        adata.obs["_group"] = "All"
        group_col = "_group"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        if metric not in adata.obs.columns:
            continue

        sns.violinplot(
            data=adata.obs,
            x=group_col,
            y=metric,
            ax=axes[idx],
            inner="quartile",
            palette="Set2"
        )

        axes[idx].set_title(f"{metric} - {title_suffix}", fontsize=12, fontweight='bold')
        axes[idx].set_xlabel(group_col, fontsize=10)
        axes[idx].set_ylabel(metric, fontsize=10)
        axes[idx].tick_params(axis='x', rotation=45)

    plt.tight_layout()

    # 保存图片
    plt.savefig(os.path.join(output_dir, f"QC_metrics_{title_suffix}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f"QC_metrics_{title_suffix}.pdf"), bbox_inches='tight')
    plt.close()

    LOGGER.info(f"已保存: {output_dir}/QC_metrics_{title_suffix}.png")
    LOGGER.info(f"已保存: {output_dir}/QC_metrics_{title_suffix}.pdf")


def generate_qc_statistics(adata_before, adata_after_filter, adata_after_doublet, outdir, group_col="sampleid"):
    """
    生成 QC 统计文件

    Parameters:
    -----------
    adata_before : AnnData
        QC 前的数据
    adata_after_filter : AnnData
        过滤后的数据
    adata_after_doublet : AnnData
        去除双细胞后的数据
    outdir : str
        输出目录
    group_col : str
        分组列名
    """
    LOGGER.info("生成 QC 统计文件...")

    os.makedirs(outdir, exist_ok=True)

    # 获取所有样本 ID
    if group_col in adata_before.obs.columns:
        sample_ids = adata_before.obs[group_col].unique()
    else:
        sample_ids = ["All"]

    # 1. 生成 cell_count_before_after_QC.xls
    cell_count_data = []
    for sample in sample_ids:
        if group_col in adata_before.obs.columns:
            before_data = adata_before.obs[adata_before.obs[group_col] == sample]
            filter_data = adata_after_filter.obs[adata_after_filter.obs[group_col] == sample]
            doublet_data = adata_after_doublet.obs[adata_after_doublet.obs[group_col] == sample]
        else:
            before_data = adata_before.obs
            filter_data = adata_after_filter.obs
            doublet_data = adata_after_doublet.obs

        row = {
            "sampleid": sample,
            "cell_count_beforeQC": len(before_data),
            "cell_count_afterQC": len(filter_data),
            "cell_count_after_rmdoublets": len(doublet_data),
            "nFeature_RNA_max_after_rmdoublets": doublet_data["nFeature_RNA"].max(),
            "nFeature_RNA_min_after_rmdoublets": doublet_data["nFeature_RNA"].min(),
            "nFeature_RNA_mean_after_rmdoublets": doublet_data["nFeature_RNA"].mean(),
            "nFeature_RNA_median_after_rmdoublets": doublet_data["nFeature_RNA"].median(),
            "nCount_RNA_max_after_rmdoublets": doublet_data["nCount_RNA"].max(),
            "nCount_RNA_min_after_rmdoublets": doublet_data["nCount_RNA"].min(),
            "nCount_RNA_mean_after_rmdoublets": doublet_data["nCount_RNA"].mean(),
            "nCount_RNA_median_after_rmdoublets": doublet_data["nCount_RNA"].median(),
            "percent.mito_max_after_rmdoublets": doublet_data["percent.mito"].max(),
            "percent.mito_min_after_rmdoublets": doublet_data["percent.mito"].min(),
            "percent.mito_mean_after_rmdoublets": doublet_data["percent.mito"].mean(),
            "percent.mito_median_after_rmdoublets": doublet_data["percent.mito"].median(),
            "percent.HB_max_after_rmdoublets": doublet_data["percent.HB"].max(),
            "percent.HB_min_after_rmdoublets": doublet_data["percent.HB"].min(),
            "percent.HB_mean_after_rmdoublets": doublet_data["percent.HB"].mean(),
            "percent.HB_median_after_rmdoublets": doublet_data["percent.HB"].median(),
        }
        cell_count_data.append(row)

    cell_count_df = pd.DataFrame(cell_count_data)
    cell_count_df.to_csv(os.path.join(outdir, "cell_count_before_after_QC.xls"), sep="\t", index=False)
    LOGGER.info("已保存 cell_count_before_after_QC.xls")

    # 2. 生成 cell_statitics_before_after_QC.xls
    cell_stat_data = []
    for sample in sample_ids:
        if group_col in adata_before.obs.columns:
            before_data = adata_before.obs[adata_before.obs[group_col] == sample]
            filter_data = adata_after_filter.obs[adata_after_filter.obs[group_col] == sample]
            doublet_data = adata_after_doublet.obs[adata_after_doublet.obs[group_col] == sample]
        else:
            before_data = adata_before.obs
            filter_data = adata_after_filter.obs
            doublet_data = adata_after_doublet.obs

        row = {
            "sampleid": sample,
            "cell_count_beforeQC": len(before_data),
            "nFeature_RNA_max_beforeQC": before_data["nFeature_RNA"].max() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_min_beforeQC": before_data["nFeature_RNA"].min() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_mean_beforeQC": before_data["nFeature_RNA"].mean() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_median_beforeQC": before_data["nFeature_RNA"].median() if "nFeature_RNA" in before_data.columns else 0,
            "nCount_RNA_max_beforeQC": before_data["nCount_RNA"].max() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_min_beforeQC": before_data["nCount_RNA"].min() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_mean_beforeQC": before_data["nCount_RNA"].mean() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_median_beforeQC": before_data["nCount_RNA"].median() if "nCount_RNA" in before_data.columns else 0,
            "percent.mito_max_beforeQC": before_data["percent.mito"].max() if "percent.mito" in before_data.columns else 0,
            "percent.mito_min_beforeQC": before_data["percent.mito"].min() if "percent.mito" in before_data.columns else 0,
            "percent.mito_mean_beforeQC": before_data["percent.mito"].mean() if "percent.mito" in before_data.columns else 0,
            "percent.mito_median_beforeQC": before_data["percent.mito"].median() if "percent.mito" in before_data.columns else 0,
            "percent.HB_max_beforeQC": before_data["percent.HB"].max() if "percent.HB" in before_data.columns else 0,
            "percent.HB_min_beforeQC": before_data["percent.HB"].min() if "percent.HB" in before_data.columns else 0,
            "percent.HB_mean_beforeQC": before_data["percent.HB"].mean() if "percent.HB" in before_data.columns else 0,
            "percent.HB_median_beforeQC": before_data["percent.HB"].median() if "percent.HB" in before_data.columns else 0,
            "cell_count_afterQC": len(filter_data),
            "nFeature_RNA_max_afterQC": filter_data["nFeature_RNA"].max() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_min_afterQC": filter_data["nFeature_RNA"].min() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_mean_afterQC": filter_data["nFeature_RNA"].mean() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_median_afterQC": filter_data["nFeature_RNA"].median() if "nFeature_RNA" in filter_data.columns else 0,
            "nCount_RNA_max_afterQC": filter_data["nCount_RNA"].max() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_min_afterQC": filter_data["nCount_RNA"].min() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_mean_afterQC": filter_data["nCount_RNA"].mean() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_median_afterQC": filter_data["nCount_RNA"].median() if "nCount_RNA" in filter_data.columns else 0,
            "percent.mito_max_afterQC": filter_data["percent.mito"].max() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_min_afterQC": filter_data["percent.mito"].min() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_mean_afterQC": filter_data["percent.mito"].mean() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_median_afterQC": filter_data["percent.mito"].median() if "percent.mito" in filter_data.columns else 0,
            "percent.HB_max_afterQC": filter_data["percent.HB"].max() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_min_afterQC": filter_data["percent.HB"].min() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_mean_afterQC": filter_data["percent.HB"].mean() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_median_afterQC": filter_data["percent.HB"].median() if "percent.HB" in filter_data.columns else 0,
            "cell_count_after_rmdoublets": len(doublet_data),
            "nFeature_RNA_max_after_rmdoublets": doublet_data["nFeature_RNA"].max() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_min_after_rmdoublets": doublet_data["nFeature_RNA"].min() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_mean_after_rmdoublets": doublet_data["nFeature_RNA"].mean() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_median_after_rmdoublets": doublet_data["nFeature_RNA"].median() if "nFeature_RNA" in doublet_data.columns else 0,
            "nCount_RNA_max_after_rmdoublets": doublet_data["nCount_RNA"].max() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_min_after_rmdoublets": doublet_data["nCount_RNA"].min() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_mean_after_rmdoublets": doublet_data["nCount_RNA"].mean() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_median_after_rmdoublets": doublet_data["nCount_RNA"].median() if "nCount_RNA" in doublet_data.columns else 0,
            "percent.mito_max_after_rmdoublets": doublet_data["percent.mito"].max() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_min_after_rmdoublets": doublet_data["percent.mito"].min() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_mean_after_rmdoublets": doublet_data["percent.mito"].mean() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_median_after_rmdoublets": doublet_data["percent.mito"].median() if "percent.mito" in doublet_data.columns else 0,
            "percent.HB_max_after_rmdoublets": doublet_data["percent.HB"].max() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_min_after_rmdoublets": doublet_data["percent.HB"].min() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_mean_after_rmdoublets": doublet_data["percent.HB"].mean() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_median_after_rmdoublets": doublet_data["percent.HB"].median() if "percent.HB" in doublet_data.columns else 0,
        }
        cell_stat_data.append(row)

    cell_stat_df = pd.DataFrame(cell_stat_data)
    cell_stat_df.to_csv(os.path.join(outdir, "cell_statitics_before_after_QC.xls"), sep="\t", index=False)
    LOGGER.info("已保存 cell_statitics_before_after_QC.xls")


def save_results(adata, outdir, adata_before=None, adata_after_filter=None, group_col="sampleid"):
    """
    保存结果
    """
    LOGGER.info(f"保存结果到 {outdir}...")

    # 创建输出目录
    os.makedirs(outdir, exist_ok=True)

    # 保存 h5ad
    adata.write_h5ad(os.path.join(outdir, "qc_result.h5ad"))
    LOGGER.info("已保存 qc_result.h5ad")

    # 保存 QC 指标
    adata.obs.to_csv(os.path.join(outdir, "qc_metrics.csv"), index=False)
    LOGGER.info("已保存 qc_metrics.csv")

    # 生成 QC 统计文件
    if adata_before is not None and adata_after_filter is not None:
        generate_qc_statistics(adata_before, adata_after_filter, adata, outdir, group_col=group_col)

        # 绘制 QC 小提琴图（过滤前）
        plot_qc_violin(adata_before, outdir, "beforeQC", group_col=group_col)

        # 绘制 QC 小提琴图（过滤后）
        plot_qc_violin(adata_after_filter, outdir, "afterQC", group_col=group_col)

    LOGGER.info("所有结果保存完成")


def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(
        description="第一步：质量控制（QC 过滤 + 双细胞去除）"
    )

    # 输入输出
    parser.add_argument("-i", "--input", required=True, help="输入文件路径 (.h5ad) 或包含样本目录的根目录")
    parser.add_argument("-o", "--outdir", required=True, help="输出目录")
    parser.add_argument("--samples", nargs="+", default=None, help="样本目录列表 (如: A1 A2 P1 P2)")

    # 参考基因组
    parser.add_argument("--refgenome", required=True, help="参考基因组目录路径")
    parser.add_argument("--mt_genelist", default="MT_genelist.gmt", help="线粒体基因列表文件名")
    parser.add_argument("--hb_genelist", default="HB_genelist.gmt", help="血红蛋白基因列表文件名")

    # QC 过滤参数
    parser.add_argument("--mingene", type=int, default=200, help="每个细胞最少基因数 (default: 200)")
    parser.add_argument("--maxgene", type=int, default=7500, help="每个细胞最多基因数 (default: 7500)")
    parser.add_argument("--minumi", type=int, default=1000, help="每个细胞最少 UMI 数 (default: 1000)")
    parser.add_argument("--maxumi", type=int, default=50000, help="每个细胞最多 UMI 数 (default: 50000)")
    parser.add_argument("--mtfilter", type=float, default=0.2, help="线粒体基因比例上限 (default: 0.2)")
    parser.add_argument("--hbfilter", type=float, default=0.05, help="血红蛋白基因比例上限 (default: 0.05)")

    # 双细胞检测
    parser.add_argument("--doublet_random_state", type=int, default=2025, help="双细胞检测随机种子 (default: 2025)")
    parser.add_argument("--skip_doublet_removal", action="store_true", help="跳过双细胞去除步骤")

    args = parser.parse_args()

    # 读取数据
    if args.samples:
        # 从多个样本目录读取数据
        LOGGER.info(f"从目录读取多个样本: {args.input}")
        adata = load_samples(args.input, args.samples)
    else:
        # 读取单个 h5ad 文件
        LOGGER.info(f"读取数据: {args.input}")
        adata = sc.read_h5ad(args.input)

    # 1. QC 指标计算
    adata = calculate_qc_metrics(
        adata,
        refgenome=args.refgenome,
        mt_genelist_file=args.mt_genelist,
        hb_genelist_file=args.hb_genelist
    )

    # 保存过滤前的数据（用于统计）
    adata_before = adata.copy()

    # 2. 过滤细胞
    adata_after_filter = filter_cells(
        adata,
        mingene=args.mingene,
        maxgene=args.maxgene,
        minumi=args.minumi,
        maxumi=args.maxumi,
        mtfilter=args.mtfilter,
        hbfilter=args.hbfilter
    )
    adata = adata_after_filter

    # 3. 去除双细胞（可选）
    if args.skip_doublet_removal:
        LOGGER.info("跳过双细胞去除步骤")
    else:
        adata = remove_doublets(adata, random_state=args.doublet_random_state)

    # 4. 保存结果（包含 QC 统计）
    save_results(
        adata,
        args.outdir,
        adata_before=adata_before,
        adata_after_filter=adata_after_filter,
        group_col="sampleid" if "sampleid" in adata.obs.columns else None
    )

    LOGGER.info("QC 分析完成！")


if __name__ == "__main__":
    main()
