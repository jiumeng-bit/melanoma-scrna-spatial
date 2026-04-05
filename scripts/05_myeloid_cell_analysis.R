#==============================================================================
# Script 5: Myeloid Cell Analysis
# Description: Subclustering, metabolic analysis, TF activity
# Corresponds to: Results 3.4
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
# 1. Myeloid Subclustering (Fig 4A, 4B)
#==============================================================================

scRNA_Mye <- subset(scRNA, celltype == 'Myeloid cell')
scRNA_Mye <- FindVariableFeatures(scRNA_Mye, selection.method = "vst", nfeatures = 2000)
scRNA_Mye <- ScaleData(scRNA_Mye)
scRNA_Mye <- RunPCA(scRNA_Mye, features = VariableFeatures(scRNA_Mye))
scRNA_Mye <- FindNeighbors(scRNA_Mye, dims = 1:10)
scRNA_Mye <- FindClusters(scRNA_Mye, resolution = 0.6)
scRNA_Mye <- RunUMAP(scRNA_Mye, dims = 1:10)

# Annotate myeloid subtypes
Mye_celltype <- c("cDC1_CLEC9A", "cDC2_CD1C", "Mono_MKI67", 
                  "Mφ_C1QB", "Mφ_CXCL8", "Mφ_CCL5")
Idents(scRNA_Mye) <- scRNA_Mye$seurat_clusters
names(Mye_celltype) <- levels(scRNA_Mye)
scRNA_Mye <- RenameIdents(scRNA_Mye, Mye_celltype)
scRNA_Mye$Mye_celltype <- Idents(scRNA_Mye)

#==============================================================================
# 2. Metabolic Analysis (Fig 4C)
#==============================================================================

scRNA_Mye <- scMetabolism::sc.metabolism(
  SeuratObj = scRNA_Mye,
  method = "VISION",
  imputation = F,
  ncores = 4,
  metabolism.type = "KEGG"
)

# Compare glycolysis across myeloid subtypes
p1 <- VlnPlot(scRNA_Mye, features = "Glycolysis / Gluconeogenesis",
              group.by = "Mye_celltype", pt.size = 0) +
  stat_compare_means(method = "anova")

#==============================================================================
# 3. Functional Markers (Fig 4D)
#==============================================================================

markers <- c("ARG1", "MRC1", "CD274", "CX3CR1", "IL10", "TGFB1")
p2 <- DotPlot(scRNA_Mye, features = markers, group.by = "Mye_celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#==============================================================================
# 4. Trajectory Analysis (Monocle3) (Fig 4E)
#==============================================================================

cds_mye <- new_cell_data_set(
  expression_data = GetAssayData(scRNA_Mye, layer = "counts"),
  cell_metadata = scRNA_Mye@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(scRNA_Mye),
                             row.names = rownames(scRNA_Mye))
)

cds_mye <- preprocess_cds(cds_mye, num_dim = 20)
cds_mye <- reduce_dimension(cds_mye, reduction_method = "UMAP")
cds_mye <- cluster_cells(cds_mye)
cds_mye <- learn_graph(cds_mye)
cds_mye <- order_cells(cds_mye)

# Plot metabolic changes along pseudotime
pseudotime_metabolism <- data.frame(
  pseudotime = pseudotime(cds_mye),
  glycolysis = scRNA_Mye@assays$RNA$count["Glycolysis / Gluconeogenesis"],
  oxidative_phosphorylation = scRNA_Mye@assays$RNA$count["Oxidative phosphorylation"],
  celltype = scRNA_Mye$Mye_celltype
)

# Plot metabolic changes
p3 <- ggplot(pseudotime_metabolism, aes(x = pseudotime)) +
  geom_smooth(aes(y = glycolysis, color = "Glycolysis"), method = "loess") +
  geom_smooth(aes(y = oxidative_phosphorylation, color = "OXPHOS"), method = "loess") +
  theme_bw() +
  labs(title = "Metabolic Changes Along Myeloid Trajectory")

#==============================================================================
# 5. pySCENIC TF Analysis
#==============================================================================

# Export data for pySCENIC (Python)
counts_export <- GetAssayData(scRNA_Mye, layer = "counts")
write.csv(counts_export, "results/myeloid_counts_for_pyscenic.csv")

cat("# Run pySCENIC:
# pyscenic grn myeloid_counts_for_pyscenic.csv tf_list.txt -o adj.csv --num_workers 8
# pyscenic ctx adj.csv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs.tbl --expression_mtx_fname melanocyte_counts_for_pyscenic.csv --mode 'dask_multiprocessing' --output reg.csv --num_workers 8
# pyscenic aucell melanocyte_counts_for_pyscenic.csv reg.csv -o aucell.csv --num_workers 8
")

saveRDS(scRNA_Mye, "results/05_myeloid_analysis.rds")
