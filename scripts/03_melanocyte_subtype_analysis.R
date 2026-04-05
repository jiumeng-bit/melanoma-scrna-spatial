#==============================================================================
# Script 3: Melanocyte Subtype Analysis
# Description: Subclustering, CNV analysis, cNMF metaprograms, metabolism
# Corresponds to: Results 3.2
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(infercnv)
  library(cNMF)
  library(scMetabolism)
  library(reticulate)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

scRNA <- readRDS("results/02_scRNA_scored.rds")

#==============================================================================
# 1. Melanocyte Subclustering (Fig 2B, 2C)
#==============================================================================

scRNA_Mel <- subset(scRNA, celltype == 'Melanocyte')
scRNA_Mel <- FindVariableFeatures(scRNA_Mel, selection.method = "vst", nfeatures = 2000)
scRNA_Mel <- ScaleData(scRNA_Mel, features = VariableFeatures(scRNA_Mel))
scRNA_Mel <- RunPCA(scRNA_Mel, features = VariableFeatures(scRNA_Mel))
scRNA_Mel <- RunHarmony(scRNA_Mel, group.by.vars = c("dataset","sample"))
scRNA_Mel <- FindNeighbors(scRNA_Mel, reduction = "harmony", dims = 1:10)
scRNA_Mel <- FindClusters(scRNA_Mel, resolution = 0.1)
scRNA_Mel <- RunUMAP(scRNA_Mel, reduction = "harmony", dims = 1:10)

# Annotate subtypes
Mel_celltype <- c("M1_SPARCL1", "M2_CRABP1", "M3_ABCB5", "M4_ISG15", 
                  "M5_ASTN2", "M6_MKI67")
Idents(scRNA_Mel) <- scRNA_Mel$seurat_clusters
names(Mel_celltype) <- levels(scRNA_Mel)
scRNA_Mel <- RenameIdents(scRNA_Mel, Mel_celltype)
scRNA_Mel$Mel_celltype <- Idents(scRNA_Mel)

#==============================================================================
# 2. CNV Analysis (inferCNV) (Fig 2D)
#==============================================================================

# Prepare annotation file
annotation_data <- data.frame(
  cell = colnames(scRNA_Mel),
  celltype = scRNA_Mel$Mel_celltype
)
write.table(annotation_data, "results/infercnv_annotation.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Create inferCNV object and run
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = GetAssayData(scRNA_Mel, layer = "counts"),
  annotations_file = "results/infercnv_annotation.txt",
  delim = "\t",
  gene_order_file = "data/gencode_v19_gene_pos.txt",
  ref_group_names = c("M2_CRABP1")
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "results/infercnv_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE
)

# Calculate CNV scores
cnv_data <- read.table("results/infercnv_output/infercnv.observations.txt", header = T, row.names = 1)
cnv_score <- colSums(cnv_data - 1)

scRNA_Mel$cnv_score <- cnv_score[colnames(scRNA_Mel)]

# Plot CNV score by subtype
p1 <- ggboxplot(scRNA_Mel@meta.data, x = "Mel_celltype", y = "cnv_score",
                fill = "Mel_celltype") +
  stat_compare_means(method = "anova")

#==============================================================================
# 3. cNMF Analysis for Metaprograms (Fig 2F, 2G, 2H)
#==============================================================================

# Export count matrix for cNMF
write.csv(GetAssayData(scRNA_Mel, layer = "counts"),
          "results/melanocyte_counts.csv")

# Run cNMF (Python code provided for reference)
cat("# Run cNMF in Python:
# import scarp as sc
# adata = sc.read_csv('results/melanocyte_counts.csv')
# sc.tl.cnmf(adata, components=7, n_iter=100)
# sc.pl.cnmf_usage(adata)
")

# Load cNMF results (after Python execution)
# cNMF_results <- read.csv("results/cnmf_usage.csv")

#==============================================================================
# 4. scMetabolism Analysis (Fig 2E, 2I)
#==============================================================================

scRNA_Mel <- scMetabolism::sc.metabolism(
  SeuratObj = scRNA_Mel,
  method = "VISION",
  imputation = F,
  ncores = 4,
  metabolism.type = "KEGG"
)

# Plot glycolysis pathway
p2 <- DimPlot.scMetabolism(scRNA_Mel, 
                           pathway = "Glycolysis / Gluconeogenesis",
                           dimention.reduction.type = "umap")

#==============================================================================
# 5. scFEA Metabolic Flux Analysis (Fig 2J)
#==============================================================================

# scFEA is Python-based - save data for external analysis
scRNA_Mel_counts <- GetAssayData(scRNA_Mel, layer = "counts")
write.table(scRNA_Mel_counts, "results/melanocyte_scfea_input.csv", 
            sep = ",", quote = FALSE)

saveRDS(scRNA_Mel, "results/03_melanocyte_subtype.rds")
