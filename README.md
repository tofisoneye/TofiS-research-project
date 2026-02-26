# Artery scRNA-seq Analysis Project

## Project Description
This project focuses on identifying cell-type-specific gene expression patterns in artery samples. The goal is to map these patterns to GWAS risk loci for temporal arteritis and other vascular diseases using the ECLIPSER tool.

## Technologies Used
* R (Seurat, dplyr)

## Repository Structure
* `First_draft.R`: Primary R script for data processing, Seurat clustering, and gene filtering.
* `README.md`: Project documentation and setup instructions.

## Setup Instructions
1. Install R (version 4.0 or higher) and RStudio.
2. Install the Seurat package: `install.packages('Seurat')`.
3. Clone this repository: `git clone https://github.com/tofisoneye/TofiS-research-project.git`.

## How to Run the Project
1. Open `First_draft.R` in RStudio.
2. Load your artery scRNA-seq RDS object (any publicly available scRNA-seq data will work).
3. Execute the script to perform normalization, clustering (clusters 0-9), and marker identification.
4. The script will output a filtered list of protein-coding genes for downstream disease mapping.
