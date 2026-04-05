#==============================================================================
# Script 1: Data Acquisition and Preprocessing
# Description: Quality control, normalization, batch correction, and cell annotation
# Corresponds to: Materials and methods 2.1
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(dplyr)
  library(future)
  library(future.apply)
  library(harmony)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(cowplot)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(SingleCellExperiment)
})

plan(multisession, workers = 8)
options(future.globals.maxSize = 20000 * 1024^2)

# Color schemes
cors <- pal_npg()(10)
cors_group <- c("#3C5488FF", "#BB0021FF")

#==============================================================================
# 1. Data Loading
#==============================================================================

cat("=== Loading datasets ===\n")

# GSE72056
scRNA1 <- data.table::fread("data/Melanoma_GSE72056.txt.gz")
clinical1 <- scRNA1[c(2,3),] %>% column_to_rownames(var = "Cell") %>% t() %>% as.data.frame()
scRNA1 <- scRNA1[-c(1,2,3),]
scRNA1 <- scRNA1[!duplicated(scRNA1$Cell),]
scRNA1 <- column_to_rownames(scRNA1, var = "Cell")
scRNA1 <- CreateSeuratObject(counts = scRNA1, project = "GSE72056", 
                             min.cells = 10, min.features = 200)
scRNA1$sample <- str_split(colnames(scRNA1), "_|-", simplify = T)[,1]
scRNA1$sample <- toupper(str_sub(scRNA1$sample, 1, 4))
scRNA1$group <- 'Melanoma'
Idents(scRNA1) <- scRNA1$sample
group_Met <- c('Primary','Primary','Metastasis','Metastasis','Metastasis',
               'Metastasis','Metastasis','Metastasis','Metastasis','Primary',
               'Metastasis','Metastasis','Metastasis','Metastasis','Metastasis',
               'Primary','Metastasis','Metastasis','Metastasis')
names(group_Met) <- levels(scRNA1)
scRNA1 <- RenameIdents(scRNA1, group_Met)
scRNA1$group_M <- Idents(scRNA1)

# GSE115978
scRNA3 <- readRDS("data/Melanoma_GSE115978.rds")
scRNA3@meta.data <- scRNA3@meta.data[,c(1:3)]
scRNA3$sample <- factor(scRNA3$orig.ident)
scRNA3$group <- 'Melanoma'
Idents(scRNA3) <- scRNA3$sample
group_Met <- c('Metastasis','Primary','Metastasis','Metastasis','Metastasis',
               'Primary','Primary','Metastasis','Metastasis','Metastasis',
               'Metastasis','Metastasis','Metastasis','Metastasis','Primary','Metastasis')
names(group_Met) <- levels(scRNA3)
scRNA3 <- RenameIdents(scRNA3, group_Met)
scRNA3$group_M <- Idents(scRNA3)

# GSE215120
scMelanoma1 <- Read10X_h5('data/Melanoma_GSE215120/GSM6622299_CM1_filtered_feature_bc_matrix.h5')
scMelanoma1 <- CreateSeuratObject(counts = scMelanoma1, project = "CM01", min.cells = 10, min.features = 200)
scMelanoma2 <- Read10X_h5('data/Melanoma_GSE215120/GSM6622300_CM2_filtered_feature_bc_matrix.h5')
scMelanoma2 <- CreateSeuratObject(counts = scMelanoma2, project = "CM02", min.cells = 10, min.features = 200)
scMelanoma3 <- Read10X_h5('data/Melanoma_GSE215120/GSM6622301_CM3_filtered_feature_bc_matrix.h5')
scMelanoma3 <- CreateSeuratObject(counts = scMelanoma3, project = "CM03", min.cells = 10, min.features = 200)
scMelanoma4 <- Read10X_h5('data/Melanoma_GSE215120/GSM6622302_CM1_lym_filtered_feature_bc_matrix.h5')
scMelanoma4 <- CreateSeuratObject(counts = scMelanoma4, project = "CM01_lym", min.cells = 10, min.features = 200)
scMelanoma <- merge(scMelanoma1, y = c(scMelanoma2, scMelanoma3, scMelanoma4))
rm(scMelanoma1, scMelanoma2, scMelanoma3, scMelanoma4)
gc()
scMelanoma$sample <- factor(scMelanoma$orig.ident)
scMelanoma$group <- 'Melanoma'
scMelanoma$group_M <- ifelse(grepl('lym', scMelanoma$sample), 'Metastasis', 'Primary')

# Healthy skin (PRJCA006797)
skin_healthy <- subset(scRNA, group == 'Healthy')
skin_healthy@meta.data <- skin_healthy@meta.data[,c(1:4)]
skin_healthy$sample <- skin_healthy$orig.ident

#==============================================================================
# 2. Merge and QC
#==============================================================================

scRNA <- merge(scRNA1, y = c(scRNA3, scMelanoma, skin_healthy))
table(scRNA$sample)
scRNA$dataset <- case_when(
  grepl('CM', scRNA$sample) ~ 'GSE215120',
  grepl('CY', scRNA$sample) ~ 'GSE72056',
  grepl('Mel', scRNA$sample) ~ 'GSE115978',
  TRUE ~ 'PRJCA006797'
)

# Quality Control
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA <- subset(scRNA, subset = nCount_RNA > 500 & nCount_RNA < 10000 &
                  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

#==============================================================================
# 3. Normalization and Batch Correction
#==============================================================================

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# Harmony batch correction
scRNA <- RunHarmony(scRNA, group.by.vars = c("sample"))
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, reduction = 'harmony', dims = 1:10)

#==============================================================================
# 4. Cell Type Annotation
#==============================================================================

library(COSG)
marker_cosg <- cosg(scRNA, groups = 'all', assay = 'RNA', slot = 'data',
                    mu = 1, remove_lowly_expressed = TRUE, expressed_pct = 0.1, n_genes_user = 100)

# Manual annotation based on markers
celltype <- c('Melanocyte', 'Keratinocyte', 'T cell', 'T cell', 'Melanocyte',
              'Melanocyte', 'Keratinocyte', 'B cell', 'Myeloid cell', 'Proliferating',
              'SMC', 'Langerhans', 'T cell', 'Fibroblast', 'Endothelial',
              'T cell', 'T cell', 'T cell', 'Melanocyte', 'Unknown', 'Proliferating')

Idents(scRNA) <- scRNA$seurat_clusters
names(celltype) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, celltype)
scRNA$celltype <- Idents(scRNA)
scRNA <- subset(scRNA, celltype != 'Unknown')

cell_order <- c("Melanocyte", "Keratinocyte", "Fibroblast", "Endothelial", "SMC",
                "T cell", "B cell", "Myeloid cell", "Langerhans")
scRNA$celltype <- factor(scRNA$celltype, levels = cell_order)

saveRDS(scRNA, "results/01_scRNA_annotated.rds")
