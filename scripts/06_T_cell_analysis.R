#==============================================================================
# Script 6: T Cell Analysis
# Description: Subclustering, CD8+ Tact analysis, exhaustion markers
# Corresponds to: Results 3.5
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(scMetabolism)
  library(Monocle3)
})

scRNA <- readRDS("results/02_scRNA_scored.rds")

#==============================================================================
# 1. T Cell Subclustering (Fig 5A, 5B)
#==============================================================================

scRNA_T <- subset(scRNA, celltype == 'T cell')
scRNA_T <- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 2000)
scRNA_T <- ScaleData(scRNA_T)
scRNA_T <- RunPCA(scRNA_T, features = VariableFeatures(scRNA_T))
scRNA_T <- FindNeighbors(scRNA_T, dims = 1:10)
scRNA_T <- FindClusters(scRNA_T, resolution = 0.8)
scRNA_T <- RunUMAP(scRNA_T, dims = 1:10)

# Annotate T cell subtypes
T_celltype <- c("CD4_Treg", "CD8_Tcyt", "CD8_Tact", "CD8_Tpro", 
                "CD8_Tifn", "CD8_Tex", "Naive_T", "NKT", "Tcm")
Idents(scRNA_T) <- scRNA_T$seurat_clusters
names(T_celltype) <- levels(scRNA_T)
scRNA_T <- RenameIdents(scRNA_T, T_celltype)
scRNA_T$T_celltype <- Idents(scRNA_T)

#==============================================================================
# 2. CD8+ Tact Glycolysis Analysis (Fig 5C)
#==============================================================================

scRNA_T <- scMetabolism::sc.metabolism(
  SeuratObj = scRNA_T,
  method = "VISION",
  imputation = F,
  ncores = 4,
  metabolism.type = "KEGG"
)

p1 <- VlnPlot(scRNA_T, features = "Glycolysis / Gluconeogenesis",
              group.by = "T_celltype", pt.size = 0) +
  stat_compare_means(method = "kruskal.test")

#==============================================================================
# 3. Co-stimulatory and Inhibitory Markers (Fig 5D, 5E, 5F)
#==============================================================================

costim_markers <- c("CD28", "TNFRSF14", "TNFRSF9")
inhibitory_markers <- c("PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4")
cytotoxic_markers <- c("GZMA", "GZMB", "PRF1", "IFNG")

p2 <- DotPlot(scRNA_T, features = c(costim_markers, inhibitory_markers),
              group.by = "T_celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#==============================================================================
# 4. Correlation Analysis (Fig 5G)
#==============================================================================

# Glycolysis vs cytotoxicity correlation
cd8_cells <- subset(scRNA_T, T_celltype %in% c("CD8_Tcyt", "CD8_Tact", "CD8_Tex"))
correlation_data <- data.frame(
  glycolysis = cd8_cells$`Glycolysis / Gluconeogenesis`,
  cytotoxicity = colMeans(GetAssayData(cd8_cells, layer = "data")[cytotoxic_markers,]),
  exhaustion = colMeans(GetAssayData(cd8_cells, layer = "data")[inhibitory_markers,])
)

cor_glyco_cytotox <- cor(correlation_data$glycolysis, correlation_data$cytotoxicity)
cor_glyco_exhaust <- cor(correlation_data$glycolysis, correlation_data$exhaustion)

#==============================================================================
# 5. Pseudotime Analysis (Fig S5)
#==============================================================================

cds_t <- new_cell_data_set(
  expression_data = GetAssayData(cd8_cells, layer = "counts"),
  cell_metadata = cd8_cells@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(cd8_cells),
                             row.names = rownames(cd8_cells))
)

cds_t <- preprocess_cds(cds_t, num_dim = 20)
cds_t <- reduce_dimension(cds_t, reduction_method = "UMAP")
cds_t <- cluster_cells(cds_t)
cds_t <- learn_graph(cds_t)
cds_t <- order_cells(cds_t)

# Track marker expression along pseudotime
pseudotime_t <- data.frame(
  pseudotime = pseudotime(cds_t),
  celltype = cd8_cells$T_celltype,
  glycolysis = cd8_cells$`Glycolysis / Gluconeogenesis`,
  IFNG = GetAssayData(cd8_cells, layer = "data")["IFNG",],
  PDCD1 = GetAssayData(cd8_cells, layer = "data")["PDCD1",],
  LAG3 = GetAssayData(cd8_cells, layer = "data")["LAG3",]
)

#==============================================================================
# 6. Chemokine Receptors (Fig 5H)
#==============================================================================

chemokine_receptors <- c("CCR5", "CXCR3", "CXCR5", "CXCR6", "CCR7")
p3 <- DotPlot(cd8_cells, features = chemokine_receptors,
              group.by = "T_celltype") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")

saveRDS(scRNA_T, "results/06_T_cell_analysis.rds")
