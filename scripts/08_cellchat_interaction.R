#==============================================================================
# Script 8: Cell-Cell Communication Analysis
# Description: CellChat analysis for M5-myeloid and Mφ_CCL5-CD8+ Tact interactions
# Corresponds to: Results 3.4, 3.5 (MIF pathway, CXCL16-CXCR6, CCL5-CCR5)
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})

#==============================================================================
# 1. Prepare Data (Fig 4F, 4G, 5I, 5J)
#==============================================================================

scRNA <- readRDS("results/02_scRNA_scored.rds")

# Subset data for cell types of interest
cell_types_interest <- c("Melanocyte", "Myeloid cell", "T cell", "B cell", 
                         "Fibroblast", "Endothelial")
scRNA_sub <- subset(scRNA, celltype %in% cell_types_interest)

# Split by condition
normaldata <- subset(scRNA_sub, group == "Healthy")
tumordata <- subset(scRNA_sub, group == "Melanoma")

#==============================================================================
# 2. Create CellChat Objects
#==============================================================================

# Healthy
normal_cellchat <- createCellChat(object = GetAssayData(normaldata, layer = "data"),
                                  meta = normaldata@meta.data,
                                  group.by = "celltype")

normal_cellchatDB <- CellChatDB.human
normal_cellchat@DB <- normal_cellchatDB
normal_cellchat <- subsetData(normal_cellchat)
normal_cellchat <- identifyOverExpressedGenes(normal_cellchat)
normal_cellchat <- identifyOverExpressedInteractions(normal_cellchat)
normal_cellchat <- projectData(normal_cellchat, PPI.human)
normal_cellchat <- computeCommunProb(normal_cellchat, raw.use = FALSE)
normal_cellchat <- filterCommunication(normal_cellchat, min.cells = 10)
normal_cellchat <- computeCommunProbPathway(normal_cellchat)
normal_cellchat <- aggregateNet(normal_cellchat)

# Melanoma
tumor_cellchat <- createCellChat(object = GetAssayData(tumordata, layer = "data"),
                                 meta = tumordata@meta.data,
                                 group.by = "celltype")

tumor_cellchatDB <- CellChatDB.human
tumor_cellchat@DB <- tumor_cellchatDB
tumor_cellchat <- subsetData(tumor_cellchat)
tumor_cellchat <- identifyOverExpressedGenes(tumor_cellchat)
tumor_cellchat <- identifyOverExpressedInteractions(tumor_cellchat)
tumor_cellchat <- projectData(tumor_cellchat, PPI.human)
tumor_cellchat <- computeCommunProb(tumor_cellchat, raw.use = FALSE)
tumor_cellchat <- filterCommunication(tumor_cellchat, min.cells = 10)
tumor_cellchat <- computeCommunProbPathway(tumor_cellchat)
tumor_cellchat <- aggregateNet(tumor_cellchat)

#==============================================================================
# 3. Compare Communication (Fig S4C, S4D)
#==============================================================================

object.list <- list(Healthy = normal_cellchat, Melanoma = tumor_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Compare interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

#==============================================================================
# 4. MIF Pathway Analysis (Fig 4F, 4G)
#==============================================================================

# MIF pathway (Melanocyte to Myeloid)
p1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
p2 <- netVisual_aggregate(cellchat, signaling = "MIF", 
                          sources.use = c("Melanocyte"),
                          targets.use = c("Myeloid cell"))

# Extract MIF pathway info
pathways.show <- c("MIF")
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, 
                            geneLR.return = FALSE)

# MIF-(CD74+CXCR4) and MIF-(CD74+CD44)
p3 <- netVisual_bubble(cellchat, sources.use = "Melanocyte",
                       targets.use = "Myeloid cell",
                       signaling = "MIF",
                       comparison = c(1, 2))

#==============================================================================
# 5. CXCL16-CXCR6 and CCL5-CCR5 (Fig 5I, 5J)
#==============================================================================

# CXCL pathway
p4 <- netVisual_bubble(cellchat, sources.use = "Myeloid cell",
                       targets.use = "T cell",
                       signaling = "CXCL",
                       comparison = c(1, 2))

# CCL pathway  
p5 <- netVisual_bubble(cellchat, sources.use = "Myeloid cell",
                       targets.use = "T cell", 
                       signaling = "CCL",
                       comparison = c(1, 2))

# Specific pairs
p6 <- netVisual_aggregate(cellchat, signaling = "CXCL",
                          sources.use = c("Myeloid cell"),
                          targets.use = c("T cell"))

#==============================================================================
# 6. M5-Myeloid Specific Analysis (using annotated subtypes)
#==============================================================================

# Load melanocyte subtypes
scRNA_Mel <- readRDS("results/03_melanocyte_subtype.rds")

# Merge annotations back
scRNA$subcelltype <- scRNA$celltype
scRNA$subcelltype[colnames(scRNA_Mel)] <- as.character(scRNA_Mel$Mel_celltype)

# Create CellChat with subtypes
scRNA_sub <- subset(scRNA, subcelltype %in% c("M5_ASTN2", "Mφ_CCL5", "Mφ_C1QB", "CD8_Tact"))

cellchat_sub <- createCellChat(object = GetAssayData(scRNA_sub, layer = "data"),
                               meta = scRNA_sub@meta.data,
                               group.by = "subcelltype")

cellchat_sub@DB <- CellChatDB.human
cellchat_sub <- subsetData(cellchat_sub)
cellchat_sub <- identifyOverExpressedGenes(cellchat_sub)
cellchat_sub <- identifyOverExpressedInteractions(cellchat_sub)
cellchat_sub <- computeCommunProb(cellchat_sub)
cellchat_sub <- computeCommunProbPathway(cellchat_sub)
cellchat_sub <- aggregateNet(cellchat_sub)

# Visualize M5 interactions
p7 <- netVisual_circle(cellchat_sub@net$weight, 
                       sources.use = "M5_ASTN2",
                       targets.use = c("Mφ_CCL5", "Mφ_C1QB"),
                       weight.scale = T)

saveRDS(cellchat, "results/08_cellchat_analysis.rds")
