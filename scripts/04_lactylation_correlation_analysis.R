#==============================================================================
# Script 4: Lactylation and Glycolysis Correlation
# Description: Pathway correlation, trajectory analysis, WGCNA
# Corresponds to: Results 3.3
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(Monocle3)
  library(dplyr)
  library(ggplot2)
  library(ggpubl)
  library(GSVA)
  library(AUCell)
  library(ComplexHeatmap)
})

scRNA_Mel <- readRDS("results/03_melanocyte_subtype.rds")

#==============================================================================
# 1. Lactylation vs Pathway Correlation (Fig 3A)
#==============================================================================

# Calculate pathway scores using GSVA
pathway_genes <- list(
  Glycolysis = c("HK2", "PFKP", "PFKL", "PGK1", "ENO1", "PKM", "LDHA"),
  mTOR = c("MTOR", "RPS6KB1", "EIF4EBP1", "RPTOR", "MLST8"),
  TGF_beta = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3"),
  EMT = c("VIM", "CDH2", "SNAI1", "TWIST1", "ZEB1")
)

expr_matrix <- as.matrix(GetAssayData(scRNA_Mel, layer = "data"))
pathway_scores <- gsva(expr_matrix, pathway_genes, method = "gsva")

# Calculate correlation with lactylation score
correlations <- sapply(rownames(pathway_scores), function(pathway) {
  cor(pathway_scores[pathway,], scRNA_Mel$lactylation_score, method = "spearman")
})

# Correlation plot
p1 <- barplot(correlations, las = 2, main = "Lactylation - Pathway Correlations")

#==============================================================================
# 2. Pseudotime Analysis (Monocle3) (Fig 3B, 3C, 3D, 3E, 3F)
#==============================================================================

# Create Monocle3 object
cds <- new_cell_data_set(
  expression_data = GetAssayData(scRNA_Mel, layer = "counts"),
  cell_metadata = scRNA_Mel@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(scRNA_Mel), 
                             row.names = rownames(scRNA_Mel))
)

# Preprocess
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Order cells by pseudotime
# Select M2 as root (healthy reference)
get_earliest_principal_node <- function(cds, cell_type = "M2_CRABP1") {
  cell_ids <- which(colData(cds)$Mel_celltype == cell_type)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(
    which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

# Plot trajectory
p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, label_branch_points = FALSE)

p3 <- plot_cells(cds, color_cells_by = "Mel_celltype",
                 label_groups_by_cluster = FALSE)

#==============================================================================
# 3. Lactylation Along Pseudotime (Fig 3E)
#==============================================================================

pseudotime_df <- data.frame(
  pseudotime = pseudotime(cds),
  lactylation_score = scRNA_Mel$lactylation_score,
  celltype = scRNA_Mel$Mel_celltype
)

p4 <- ggplot(pseudotime_df, aes(x = pseudotime, y = lactylation_score)) +
  geom_smooth(method = "loess", color = "red") +
  geom_point(aes(color = celltype), alpha = 0.5) +
  theme_bw() +
  labs(title = "Lactylation Score Along Pseudotime")

#==============================================================================
# 4. hdWGCNA Analysis (Fig 3G, 3H)
#==============================================================================

# Note: hdWGCNA is a separate R package
# Install from GitHub: devtools::install_github('smorabit/hdWGCNA')

library(hdWGCNA)

# Set up WGCNA
scRNA_Mel <- SetupForWGCNA(
  scRNA_Mel,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Melanocyte"
)

# Construct metacells
scRNA_Mel <- MetacellsByGroups(
  seurat_obj = scRNA_Mel,
  group.by = c("Mel_celltype", "sample"),
  k = 25,
  max_shared = 10
)

# Normalize metacells
scRNA_Mel <- NormalizeMetacells(scRNA_Mel)

# Set soft threshold
scRNA_Mel <- SetDatExpr(
  scRNA_Mel,
  group_name = "Mel_celltype",
  group.by = "Mel_celltype"
)

# Test soft powers
soft_power <- TestSoftPowers(scRNA_Mel, networkType = "signed")

# Construct network
scRNA_Mel <- ConstructNetwork(scRNA_Mel, soft_power = 5)

# Compute modules
scRNA_Mel <- ModuleEigengenes(scRNA_Mel)
scRNA_Mel <- ModuleConnectivity(scRNA_Mel)

# Correlate with lactylation score
scRNA_Mel <- ModuleTraitCorrelation(
  scRNA_Mel,
  traits = "lactylation_score"
)

# Get lactylation-related modules (cor > 0.4, p < 0.05)
module_trait <- GetModuleTraitCorrelation(scRNA_Mel)
lactylation_modules <- names(which(module_trait$cor[,"lactylation_score"] > 0.4 &
                                     module_trait$pval[,"lactylation_score"] < 0.05))

saveRDS(scRNA_Mel, "results/04_melanocyte_lactylation.rds")
