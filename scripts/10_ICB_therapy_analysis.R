#==============================================================================
# Script 10: ICB Therapy Response Analysis
# Description: Analysis of immunotherapy response in GSE120575 and IMvigor210
# Corresponds to: Results 3.8 (Fig 9)
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(IMvigor210CoreBiologies)
  library(dplyr)
  library(pheatmap)
})

#==============================================================================
# 1. GSE120575 scRNA-seq Analysis (Fig 9A, 9B, 9C, 9D, 9E)
#==============================================================================

# Load GSE120575 data (ICB-treated melanoma)
icb_data <- readRDS("data/GSE120575_seurat.rds")

# Annotate cell types (already annotated in original study)
icb_data$cell_type <- icb_data$cellType

# Calculate lactylation score
lactylation_genes <- c("HIF1A", "LDHA", "ENO1", "GAPDH", "PGK1", "PGAM1", "PKM", "HK2")
icb_data <- AddModuleScore(icb_data, features = list(lactylation_genes),
                           name = "lactylation_score")

# Cell type proportions by response
p1 <- table(icb_data$cell_type, icb_data$response) %>%
  as.data.frame() %>%
  ggplot(aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Response", y = "Proportion", fill = "Cell Type")

# Lactylation scores by response (baseline and post-treatment)
baseline_data <- subset(icb_data, timepoint == "Pre")
post_data <- subset(icb_data, timepoint == "On")

# Baseline comparison
p2 <- ggboxplot(baseline_data@meta.data, x = "response", y = "lactylation_score1",
                fill = "response") +
  stat_compare_means(method = "wilcox.test") +
  facet_wrap(~cell_type) +
  labs(title = "Baseline Lactylation Scores")

# SLC25A39 expression change (Fig 9E)
sc_icb <- AddModuleScore(sc_icb, features = list("SLC25A39"), name = "SLC25A39")

p3 <- sc_icb@meta.data %>%
  group_by(response, timepoint) %>%
  summarise(mean_SLC25A39 = mean(SLC25A391)) %>%
  ggplot(aes(x = timepoint, y = mean_SLC25A39, group = response, color = response)) +
  geom_line() +
  geom_point() +
  labs(title = "SLC25A39 Expression Dynamics")

#==============================================================================
# 2. IMvigor210 Cohort Analysis (Fig 9F, 9G, 9H)
#==============================================================================

# Load IMvigor210 data
data("IMvigor210CoreBiologies")

# Get expression and clinical data
cds <- IMvigor210CoreBiologies::count_mtx
clin <- IMvigor210CoreBiologies::clinical_annotation

# Calculate LRGS (using coefficients from script 09)
LRGS_model <- readRDS("results/09_LRGS_model.rds")
lrgs_genes <- LRGS_model$LRGS_genes
lrgs_coef <- LRGS_model$coefficients

# Match genes
common_genes <- intersect(lrgs_genes, rownames(cds))
cds_sub <- cds[common_genes,]

# Calculate LRGS
lrgs_score <- apply(cds_sub, 2, function(x) sum(x * lrgs_coef[common_genes]))
clin$LRGS <- lrgs_score[rownames(clin)]

# Pathway scores
glycolysis_genes <- c("HK2", "PFKP", "PGK1", "ENO1", "PKM", "LDHA")
lactylation_genes <- c("HIF1A", "LDHA", "ENO1", "GAPDH", "PGK1", "PGAM1", "PKM", "HK2")

clin$glycolysis_score <- colMeans(cds[glycolysis_genes,])
clin$lactylation_score <- colMeans(cds[lactylation_genes,])

# Compare by response (Fig 9F, 9G)
p4 <- ggboxplot(clin, x = "Best.Confirmed.Overall.Response", 
                y = "lactylation_score",
                fill = "Best.Confirmed.Overall.Response") +
  stat_compare_means(method = "kruskal.test")

p5 <- ggboxplot(clin, x = "Best.Confirmed.Overall.Response",
                y = "glycolysis_score",
                fill = "Best.Confirmed.Overall.Response") +
  stat_compare_means(method = "kruskal.test")

# Survival analysis by LRGS (Fig 9H)
res.cut <- surv_cutpoint(clin, time = "os", event = "censOS",
                         variables = "LRGS")
res.cat <- surv_categorize(res.cut)

fit <- survfit(Surv(os, censOS) ~ LRGS, data = res.cat)
p6 <- ggsurvplot(fit, data = res.cat, pval = TRUE, risk.table = TRUE,
                 title = "OS by LRGS (IMvigor210)")

# Response rate by LRGS (Supplementary Fig S9)
response_lrgs <- table(clin$LRGS, clin$Best.Confirmed.Overall.Response)
chi_test <- chisq.test(response_lrgs)

saveRDS(list(
  gse120575_results = icb_data,
  imvigor210_results = clin,
  response_analysis = response_lrgs
), "results/10_ICB_analysis.rds")
