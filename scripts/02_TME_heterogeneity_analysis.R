#==============================================================================
# Script 2: TME Heterogeneity Analysis
# Description: Cell proportion analysis, hypoxia/lactylation scoring, metabolic assessment
# Corresponds to: Results 3.1
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(AUCell)
  library(GSVA)
  library(GSEABase)
})

scRNA <- readRDS("results/01_scRNA_annotated.rds")

#==============================================================================
# 1. Cell Proportion Analysis (Fig 1D, 1E, 1F)
#==============================================================================

# Calculate cell proportions
prop_data <- scRNA@meta.data %>%
  group_by(sample, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# Proportion plot
p1 <- ggplot(prop_data, aes(x = sample, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cors) +
  theme_bw() +
  labs(title = "Cell Type Proportions")

# K-means clustering of samples by cell proportions
prop_matrix <- prop_data %>%
  select(sample, celltype, prop) %>%
  tidyr::spread(key = celltype, value = prop, fill = 0)
rownames(prop_matrix) <- prop_matrix$sample
prop_matrix$sample <- NULL

set.seed(123)
kmeans_result <- kmeans(prop_matrix, centers = 4)

#==============================================================================
# 2. Hypoxia and Lactylation Scoring (Fig 1G, 1H, 1I)
#==============================================================================

# Gene sets
hypoxia_genes <- c("HIF1A", "EPAS1", "VEGFA", "PGK1", "ENO1", "SLC2A1", "CA9")
lactylation_genes <- c("HIF1A", "LDHA", "ENO1", "GAPDH", "PGK1", "PGAM1", "PKM", "HK2")

# AUCell scoring
expr_matrix <- GetAssayData(scRNA, slot = "data")
geneSets <- GeneSet(lactylation_genes, setName = "Lactylation_Signature")
aucellRankings <- AUCell_buildRankings(exprMatrix = expr_matrix)
aucellResults <- AUCell_calcAUC(geneSets, aucellRankings, aucMaxRank = nrow(aucellRankings)*0.05)
scRNA$lactylation_score <- as.vector(getAUC(aucellResults))

# Hypoxia score
geneSets_hypoxia <- GeneSet(hypoxia_genes, setName = "Hypoxia_Signature")
aucellResults_hypoxia <- AUCell_calcAUC(geneSets_hypoxia, aucellRankings, aucMaxRank = nrow(aucellRankings)*0.05)
scRNA$hypoxia_score <- as.vector(getAUC(aucellResults_hypoxia))

# Lactylation by condition (Fig 1G)
p2 <- ggboxplot(scRNA@meta.data, x = "group", y = "lactylation_score",
                fill = "group", palette = cors_group) +
  stat_compare_means(method = "wilcox.test") +
  facet_wrap(~celltype, scales = "free")

# Hypoxia by condition (Fig 1H)
p3 <- ggboxplot(scRNA@meta.data, x = "group", y = "hypoxia_score",
                fill = "group", palette = cors_group) +
  stat_compare_means(method = "wilcox.test")

# HIF1A+ melanocytes lactylation (Fig 1I)
scRNA$HIF1A_exp <- scRNA@assays$RNA$counts['HIF1A',]
scRNA$HIF1A_pos <- ifelse(scRNA$HIF1A_exp > 0, "HIF1A+", "HIF1A-")

melanocytes <- subset(scRNA, celltype == "Melanocyte")
p4 <- ggboxplot(melanocytes@meta.data, x = "HIF1A_pos", y = "lactylation_score",
                fill = "HIF1A_pos") +
  stat_compare_means(method = "wilcox.test")

#==============================================================================
# 3. Glucose Metabolism Pathway Analysis (Fig 1J)
#==============================================================================

# GSVA for metabolic pathways
kegg_metabolism <- readRDS("data/kegg_metabolic_pathways.rds")
glycolysis_genes <- kegg_metabolism[["Glycolysis / Gluconeogenesis"]]

# GSVA analysis
expr_matrix_gsva <- as.matrix(GetAssayData(scRNA, layer = "data"))
gsva_scores <- gsva(expr_matrix_gsva, 
                    gset.idx.list = list(Glycolysis = glycolysis_genes),
                    method = "gsva", kcdf = "Gaussian")

scRNA$glycolysis_score <- gsva_scores[1,]

# Compare glycolysis between healthy and tumor
p5 <- ggboxplot(scRNA@meta.data, x = "group", y = "glycolysis_score",
                fill = "group", palette = cors_group) +
  stat_compare_means(method = "wilcox.test") +
  facet_wrap(~celltype)

saveRDS(scRNA, "results/02_scRNA_scored.rds")
