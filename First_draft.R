###Name: Tofi Soneye
#Affiliation: Segre Lab, Mass Eye and Ear, HMS
#Purpose: Figure generation for research paper

#----
#Load required libraries, packages and the dataset to be analyzed
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)

install.packages("readr")
#dittoseq for dot plots
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")

#Step 1: Load the batch-corrected TA dataset
load("~/Downloads/IntegratedData (3).RData")
#-----
##The dataset, integratedSeurat had the following preliminary work done on it:
#QC (Gupta lab previously removed cells with >5% mtDNA, <25 NCount_RNA, <3 non-zero genes)
#PCA, Scaling, UMAP,Integration (batch correction),Clustering

#My next steps: More stringent  QC(removing cells with less than 500 genes and more than the 99.7th percentile nCount),  Differential gene expression,GSEA and cell annotations.
#in addition to visualizing and generating plots that "describe" the data

#Set the default assay to RNA (RNA assay has raw counts/ full gene list)
DefaultAssay(integratedSeurat) <- "RNA"
nFeature_cutoff <- 500
# Determine the value of the 99.7th percentile nCount_RNA
nCount_cuttof <- quantile(integratedSeurat$nCount_RNA, probs = 0.997) #nCount cutoff= 44445.75

#Step 2: QC
ta_filtered <- subset(integratedSeurat, nFeature_RNA >nFeature_cutoff & nCount_RNA < nCount_cuttof )

#Figure 1
#####
#1a. UMAP colored by disease/control

#first rename disease and control status: raw data lists them as "Neg", "Pos"
# Changing "Neg" to "Control", and "Pos" to disease

ta_filtered@meta.data$Disease[ ta_filtered@meta.data$Disease == "Neg" ] <- "Control"
ta_filtered@meta.data$Disease[ ta_filtered@meta.data$Disease == "Pos" ] <- "Disease"

#UMAP colorcoded by disease status
DimPlot(ta_filtered, reduction = "umap", group.by = "Disease") +
  ggtitle("Temporal artery data color-coded by disease status")

#UMAP colorcoded by donor ID
DimPlot(ta_filtered, reduction = "umap", group.by = "SampleID") +
  ggtitle("Temporal artery data color-coded by donor ID")

#UMAP colorcoded by cluster number
DimPlot(ta_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Temporal artery data color-coded by unique cell clusters")


#Supplementary figure 1 : UMAP of each donor vs all others ( a clearer way to observe the distributuon between all the donors that make up the data)
######
#AT16
ta_filtered$label <- ifelse(ta_filtered$SampleID== "AT16", "AT16","others")
AT16 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor AT16")

#AT11
ta_filtered$label <- ifelse(ta_filtered$SampleID== "AT11", "AT11","others")
AT11 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor AT11")

#AT18
ta_filtered$label <- ifelse(ta_filtered$SampleID== "AT18", "AT18","others")
AT18 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor AT18")

#DRGA001
ta_filtered$label <- ifelse(ta_filtered$SampleID== "DRGA001", "DRGA001","others")
DRGA001 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor DRGA001")

#DRGA002
ta_filtered$label <- ifelse(ta_filtered$SampleID== "DRGA002", "DRGA002","others")
DRGA002 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor DRGA002")

#DRGA005
ta_filtered$label <- ifelse(ta_filtered$SampleID== "DRGA005", "DRGA005","others")
DRGA005 <-DimPlot(ta_filtered, group.by = "label", cols = c("red", "grey")) +
  ggtitle("Distribution of cells from donor DRGA005")


##FIGURE 2: BUBBLE MAP OF AVERAGE EXPRESSION OF TOP 5 GENE MARKERS PER CELL CLUSTER (All samples, disease only, control only)




#####
#cluster 0 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster0,
  vars = c("TMC8", "IKZF3", "ZAP70","RASAL3"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 0")
cluster_0_pc_genes

#cluster 1 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster1,
  vars = c("ABCA10", "VIT", "AOX1","ZNF385B","SCARA5"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 1")
cluster_1_pc_genes

#cluster 2 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster2,
  vars = c("TBX2", "ITGA11", "KRT7", "TRPC4", "ANO1"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 2")
cluster_2_pc_genes

#cluster 3 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster3,
  vars = c( "CSF2RA", "MS4A7", "DCSTAMP", "ALOX15B", "TM4SF19"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 3")
cluster_3_pc_genes

#cluster 4 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster4,
  vars = c( "NOSTRIN", "TPO", "BTNL9","NR5A2","LNX1"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 4")
cluster_4_pc_genes

#cluster 5 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster5,
  vars = c( "BMP6", "CPAMD8", "CRTAC1", "SULT1B1", "EDN1"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 5")
cluster_5_pc_genes

#cluster 6 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster6,
  vars = c( "RYR2", "HPSE2", "UNC13C","COL4A6","WNK2"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 6")
cluster_6_pc_genes

#cluster 7 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster7,
  vars = c("BCL11B", "CD96", "ITK", "IL7R","GRAP2"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 7")
cluster_7_pc_genes


#cluster 8 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster8,
  vars = c( "PKHD1L1", "DLGAP2", "PIK3C2G", "TFF3", "CCL21"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 8")
cluster_8_pc_genes


#cluster 9 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster9,
  vars = c(  "NRXN1", "XKR4", "FRMD5", "SLC35F1", "NKAIN3"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 9")
cluster_9_pc_genes

#cluster 10 bubble plot
dittoSeq::dittoDotPlot(
  object = ta_cluster10,
  vars = c("ESM1", "OTC", "FAM43A", "CDH8", "ZNF365"),
  group.by = "seurat_clusters"
) +
  ggtitle("Top 5 protein-coding genes in cluster 10")
cluster_10_pc_genes




#####
##FIGURE 2: BUBBLE MAP OF AVERAGE EXPRESSION OF TOP 5 GENE MARKERS PER CELL CLUSTER (All samples, disease only, control only)
#a: all samples
#####
# Run Wilcoxon Rank Sum test
# only.pos = TRUE ensures we only get genes enriched in the cell type
all_markers <- FindAllMarkers(ta_filtered,
                              only.pos = TRUE,
                              min.pct = 0.1,
                              logfc.threshold = 0.58,
                              test.use = "wilcox")

# 2. Filter for your specific significance and top 5 genes
top20_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%      # FDR cutoff
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

# Extract gene names for the plot
genes_to_plot <- unique(top20_markers$gene)

#the gene list was taken from genes_to_plot. I manually selected the top protein coding genes
genes_to_plot <- c(
  "ZAP70", "TBC1D10C", "LYZ", "TMC8", "IKZF3",
  "ZNF385B", "SCARA5", "ANKFN1", "VIT", "KCND2",
  "KRT7", "PLA2G5", "TRPC4", "CDH2", "ACAN",
  "HTRA4", "F13A1", "MS4A7", "CLEC4E", "CSF2RA",
  "BTNL9", "TPO", "NOSTRIN", "RAB3C", "STC1",
  "CLDN10", "CPAMD8", "BMP6", "EDN1", "CRTAC1",
  "MYO3A", "MAPK4", "UNC13C", "PTCHD1", "WNK2",
  "CD28", "BCL11B", "SPOCK2", "TRBC1", "ITK",
  "TSPEAR", "PIK3C2G", "PROX1", "PKHD1L1", "CDH22",
  "SLC9A3", "GRIK3", "XKR4", "SORCS3", "SYT10",
  "ESM1", "DCTPP1", "CA8", "ZNF365", "OXSM"
)

DotPlot(ta_filtered, features = genes_to_plot, dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("lightgrey", "blue", "darkblue")) +

  labs(title = "Top 5 Protein-Coding Markers in each cluster (all samples)",
       subtitle = "Selection criteria: Adj. p-value < 0.05 & log2FC > 0.58 (FC > 1.5)",
       x = "Protein-coding Genes",  # Labels X-axis
       y = "Cell  clusters") +
  theme(axis.text.x = element_text(size = 6))



#Supplementary figure 2
#A: Bubble map of average expression of top 5 gene markers (protein coding) per cell cluster in disease samples.
#subset the filtered object based on disease status
disease_group <- subset(ta_filtered, subset = Disease == "Disease")

# only.pos = TRUE ensures we only get genes enriched in the cell type
disease_markers <- FindAllMarkers(disease_group,
                                  only.pos = TRUE,
                                  min.pct = 0.1,
                                  logfc.threshold = 0.58,
                                  test.use = "wilcox")

# 2. Filter significance and top 5 genes
top_disease_markers <- disease_markers  %>%
  filter(p_val_adj < 0.05) %>%      # FDR cutoff
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

disease_genes <- c(
  "ZAP70", "TBC1D10C", "TMC8", "IKZF3", "CCDC88B",
  "KCND2", "COL11A1", "ABCA10", "ZNF385B", "SCARA5",
  "KRT7", "TRPC4", "PLA2G5", "ACAN", "CDH2",
  "HTRA4", "F13A1", "MS4A7", "CLEC4E", "CSF2RA",
  "BTNL9", "TPO", "NOSTRIN", "RAB3C", "STC1",
  "CPAMD8", "CLDN10", "BMP6", "EDN1", "CRTAC1",
  "MYO3A", "GPR20", "UNC13C", "MAPK4", "SLC6A17",
  "CD28", "BCL11B", "TRBC1", "SPOCK2", "ICOS",
  "PKHD1L1", "PROX1", "CDH22", "TFF3", "TSPEAR",
  "SLC9A3", "PRIMA1", "GRIK3", "SYT10", "PCSK2",
  "ESM1", "CA8", "DCTPP1", "OXSM", "ACSM3"
)

#bubble plot with the gene list
DotPlot(disease_group, features = disease_genes, dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("lightgrey", "blue", "darkblue")) +

  labs(title = "Top 5 Protein-Coding Markers in each cluster (Disease group)",
       subtitle = "Selection criteria: Adj. p-value < 0.05 & log2FC > 0.58 (FC > 1.5)",
       x = "Protein-coding Genes",  # Labels X-axis
       y = "Cell  clusters") +
  theme(axis.text.x = element_text(size = 6))






#B:Bubble map of average expression of top 5 gene markers (protein coding) per cell cluster in control samples.
control_group <- subset(ta_filtered, subset = Disease == "Control")

# only.pos = TRUE ensures we only get genes enriched in the cell type
healthy_markers <- FindAllMarkers(control_group,
                                  only.pos = TRUE,
                                  min.pct = 0.1,
                                  logfc.threshold = 0.58,
                                  test.use = "wilcox")

# 2. Filter specific significance and top 5 genes
top_healthy_markers <- healthy_markers %>%
  filter(p_val_adj < 0.05) %>%      # FDR cutoff
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()


control_genes <- c(
  "SNX22", "MS4A1", "CD48", "NIBAN3", "SELL",
  "ANKFN1", "LMX1A", "CDH20", "MYOC", "ZNF385B",
  "RBM20", "INHBA", "TRPC4", "ST6GAL2", "DCDC2C",
  "MS4A6A", "FGD2", "CLEC7A", "RGS18", "IL1R2",
  "TPO", "BTNL9", "MEOX1", "NOSTRIN", "AIM2",
  "ADGRF4", "F5", "MYOZ3", "NWD1", "GDF7",
  "PTCHD1", "MYO3A", "MAPK4", "ACTG2", "AVPR1A",
  "UBASH3A", "CD28", "SLAMF1", "AGAP2", "IL7R",
  "TSPEAR", "PIK3C2G", "POF1B", "PROX1", "SPHKAP",
  "BCAS1", "SORCS3", "MPZ", "CTXND1", "GFRA3"
)

#bubble plot with the gene list
DotPlot(control_group, features = control_genes, dot.scale = 6) +
  RotatedAxis() +
  scale_colour_gradientn(colours = c("lightgrey", "blue", "darkblue")) +

  labs(title = "Top 5 Protein-Coding Markers in each cluster (Control group)",
       subtitle = "Selection criteria: Adj. p-value < 0.05 & log2FC > 0.58 (FC > 1.5)",
       x = "Protein-coding Genes",  # Labels X-axis
       y = "Cell  clusters") +
  theme(axis.text.x = element_text(size = 6))


#Tables
#1:QC of TA snRNA-seq data (gene and cell filtering).


#the tables below have markers selected by the Wilcoxon rank sum test with these parameters:
#only.pos = TRUE #only shows genes enriched in the target cluster realtive to other cells
#min.pct = 0.1 #filters for genes that are expressed in at least 10% of cells
#logfc.threshold = 0.58  #minimum LFC of selected markers corresponds to a 1.5 fold difference in expression



#2:DGE results for each cluster in all samples.
write.csv(all_markers, "DGE results for each cluster in all samples.csv")

#3:DGE results for each cluster in control samples.
write.csv(control_markers, "DGE results for each cluster in controls.csv")

#4:DGE results for each cluster in disease samples.
write.csv(disease_markers, "DGE results for each cluster in TA cells.csv")

#5:Top 20 highly expressed genes for each cluster in all samples (FDR<0.05 and highest log(fold-change))
write.csv(top20_markers, "Top 20 highly expressed genes for each cluster in all samples.csv")



####GSEA genes
#G-profiler on DEG with FDR<0.05
#1: All samples
all_gsea_markers <- FindAllMarkers(ta_filtered,
                                   only.pos = TRUE,
                                   min.pct = 0.1,
                                   logfc.threshold = 0.58, # Approx 1.5 fold change
                                   test.use = "wilcox")

# 2. Filter for FDR (p_val_adj) < 0.05
all_samples_gsea_markers <- all_gsea_markers %>%
  filter(p_val_adj < 0.05)

#save as csv
write.csv(all_samples_gsea_markers, "GSEA markers in all samples.csv")

#2: Control samples
control_markers <- FindAllMarkers(healthy_group,
                                  only.pos = TRUE,
                                  min.pct = 0.1,
                                  logfc.threshold = 0.58, # Approx 1.5 fold change
                                  test.use = "wilcox")

# 2. Filter for FDR (p_val_adj) < 0.05
sig_control_markers <- control_markers %>%
  filter(p_val_adj < 0.05)
#save as csv
write.csv(sig_control_markers, "GSEA markers in healthy samples.csv")



#3: Disease samples
Disease_markers <- FindAllMarkers(disease_group,
                                  only.pos = TRUE,
                                  min.pct = 0.1,
                                  logfc.threshold = 0.58, # Approx 1.5 fold change
                                  test.use = "wilcox")

# 2. Filter for FDR (p_val_adj) < 0.05
sig_disease_markers <- Disease_markers %>%
  filter(p_val_adj < 0.05)

write.csv(Disease_markers, "GSEA markers in disease samples.csv")

















