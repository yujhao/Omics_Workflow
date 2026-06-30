
# Omics\_Workflow

Integrated transcriptomics analysis workflows for single\-cell, bulk and spatial RNA\-seq data\.

## Repository Structure

### Single\_cell\_RNA\-seq

Complete single\-cell RNA\-seq analytical scripts under `scRNA-seq/`: Quality control, cell annotation, RNA velocity, trajectory inference, cell communication and pathway enrichment visualization\.

### Bulk\_RNA\-seq

- `bulkRNA-seq`: Standard bulk transcriptome pipeline, differential expression \& functional enrichment\.
- `ChIP-seq`: Reserved workflow for epigenetic peak analysis\.
- `common`: Shared functions for bulk modules\.

### Spatial\_transcriptomics\_RNA\-seq

- `10×Visium`: Classic spatial transcriptome workflow for spatial clustering, tissue marker visualization and spatial ligand\-receptor analysis\.
- `Stereo-seq`: Reserved workflow for Stereo\-seq spatial transcriptomics analytical pipeline.

### utils

Universal plotting, statistical tools and gene annotation resources shared by all three modules\.

## Features

- Modular folder structure, easy to expand new omics workflows
- Standardized plotting code for publication figures
- Reusable shared scripts to avoid redundant code

## Author \&amp; Acknowledgement

Maintainer: yujhao
If this workflow is adopted for academic research and publication, please cite this repository. Bioinformatics analysis was assisted by OE Biotech Co., Ltd. (Shanghai, China).
GitHub: https://github\.com/yujhao/Omics\_Workflow\.git

## License

MIT License

