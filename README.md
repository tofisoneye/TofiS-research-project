# Single-Cell RNA Sequencing Analysis for Giant Cell Arteritis (GCA) Drug Target Discovery

## Overview

This project applies machine learning and multi-omic integration to identify **cell-type-specific gene expression patterns** in temporal artery tissue samples from Giant Cell Arteritis (GCA) patients and healthy controls. By mapping disease-associated GWAS loci to these cell-type-specific signatures using the **ECLIPSR** tool, this work aims to pinpoint which cellular populations and genes drive GCA pathogenesis, enabling targeted drug therapy development.

## Project Goals

1. **Characterize cell types** in temporal artery tissue via single-cell RNA sequencing (scRNA-seq) clustering
2. **Identify disease-associated genes** through differential expression analysis (GCA vs. healthy control)
3. **Map GWAS loci to cell types** using ECLIPSR to reveal which cell populations harbor genetic drivers of GCA
4. **Generate precision medicine insights** to support targeted therapeutic development

## Background: Why This Matters

Giant Cell Arteritis (GCA) is the most common vasculitis in older adults, causing inflammation of medium and large arteries. Current understanding of cellular and genetic drivers is limited. By integrating single-cell transcriptomics with genome-wide association studies (GWAS), this project identifies which immune and vascular cell types are most affected, which genes are critical in disease pathogenesis, and which genetic variants act specifically in disease-relevant cell types. This enables precision medicine: targeting therapies to the cells and pathways that actually cause disease.

## Technologies & Tools

- **scRNA-seq Analysis:** R, Seurat (≥4.0)
- **Data Processing:** dplyr, tidyr
- **Visualization:** ggplot2
- **GWAS Integration:** ECLIPSR
- **Computing Environment:** R on HPC Terminal (Linux/Bash)

## Repository Structure


scripts/              → R code for data processing, clustering, annotation
data/                → Cell annotations, filtering results, cluster markers
outputs/             → Differential expression results, visualizations
README.md            → Project documentation
First_draft.R        → Primary Seurat analysis script


## Setup Instructions

### Prerequisites

- R (≥4.0)
- RStudio (optional but recommended)
- Linux/HPC environment for ECLIPSR pipeline

### Step 1: Install R Packages

```r
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "Matrix"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Step 2: Clone Repository

```bash
git clone https://github.com/tofisoneye/TofiS-research-project.git
cd TofiS-research-project
```

### Step 3: Prepare Data

Place temporal artery scRNA-seq data (Seurat RDS object) in the project directory or update file paths in scripts.

### Step 4: Run Analysis

```bash
# Run main analysis script
Rscript First_draft.R
```

## How to Run the Project

**Input:** Seurat RDS object with scRNA-seq data and clusters

**Execution:**
1. Open `First_draft.R` in RStudio or terminal
2. Load your temporal artery scRNA-seq data
3. Execute script to perform normalization, clustering, and marker gene identification
4. Outputs include cluster visualizations (PNGs) and differential expression results (CSVs)

**Output files:**
- Cluster plots and quality control PDFs
- CSV files with differential expression results per cluster
- Cell type annotations based on GSEA enrichment

## Cell Type Annotations

11 clusters identified via scRNA-seq:

- **Cluster 0, 3:** NKT/Macrophage (immune cells)
- **Cluster 1:** Fibroblast A/B
- **Cluster 2, 6:** Fibromyocyte/VSMC (vascular smooth muscle)
- **Cluster 4, 5, 8:** Endothelial cells (EC)
- **Cluster 7:** Natural Killer T cells (NKT)
- **Cluster 9:** Neuronal
- **Cluster 10:** Macrophage/EC (mixed)

## Differential Expression Analysis

Three DEG tables generated:
- **all_markers:** DEGs across all cells
- **disease_markers:** DEGs in GCA samples
- **control_markers:** DEGs in healthy controls

Formatted for ECLIPSR integration with standardized columns: p_val, log2FC, pct.1, pct.2, pvals_fdr, celltype, gene, tissue

## ECLIPSR Integration

DEG tables are processed and formatted for ECLIPSR pipeline to map GWAS loci to cell-type-specific gene expression, identifying which cell types harbor genetic drivers of GCA.

## Development Progress

**Completed:**
- scRNA-seq quality control and clustering (11 clusters)
- GSEA-based cell type annotation
- Differential expression analysis
- ECLIPSR input formatting

**In Progress:**
- ECLIPSR GWAS-to-cell-type mapping
- Validation against literature and public datasets

**Next Steps:**
- Manuscript preparation and writing
- Publication

## Methods

**scRNA-seq pipeline:** Quality control, normalization, PCA, clustering, UMAP visualization using Seurat

**Cell type annotation:** Marker gene identification and gene set enrichment analysis (GSEA)

**Differential expression:** Seurat FindMarkers() function comparing GCA vs. control samples per cell type

**GWAS integration:** ECLIPSR software maps GWAS summary statistics to cell-type-specific gene expression

## Data Availability

Temporal artery scRNA-seq data from GCA patients and healthy controls. Contact project lead for data access information.

## References

Analysis uses Seurat (R package for scRNA-seq), ECLIPSR (for GWAS integration), and GSEA (for cell type validation).

## Author

Tofi Soneye
