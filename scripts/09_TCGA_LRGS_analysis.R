#==============================================================================
# Script 9: TCGA Analysis and LRGS Construction
# Description: Bulk RNA analysis, survival analysis, LASSO regression
# Corresponds to: Results 3.7 (Fig 7, S8)
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(limma)
  library(sva)
  library(survival)
  library(survminer)
  library(glmnet)
  library(timeROC)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

#==============================================================================
# 1. Data Acquisition
#==============================================================================

# TCGA-SKCM
query_tcga <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# GDCdownload(query_tcga)
# tcga_data <- GDCprepare(query_tcga)

# GTEx data for normal controls
# gtex_skin <- readRDS("data/GTEx_skin.rds")

# GSE65904 (validation)
# gse65904 <- getGEO("GSE65904")

#==============================================================================
# 2. Differential Expression Analysis (Fig 7A, 7B)
#==============================================================================

# Merge TCGA and GTEx
# counts_combined <- cbind(assay(tcga_data), assay(gtex_skin))
# coldata <- data.frame(
#   condition = c(rep("Tumor", ncol(tcga_data)), rep("Normal", ncol(gtex_skin))),
#   row.names = colnames(counts_combined)
# )

# Combat batch correction
# batch <- c(rep("TCGA", ncol(tcga_data)), rep("GTEx", ncol(gtex_skin)))
# adjusted_counts <- ComBat(dat = counts_combined, batch = batch)

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = rmb_cmb,  # batch-corrected counts
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

#==============================================================================
# 3. Lactylation-related Genes from hdWGCNA
#==============================================================================

# Load WGCNA results from single-cell analysis
scRNA_Mel <- readRDS("results/04_melanocyte_lactylation.rds")
module_eigengenes <- GetModuleEigengenes(scRNA_Mel)

# Get lactylation-related modules (blue, pink, black, purple, greenyellow)
lactylation_genes <- GetWGCNAGenes(scRNA_Mel, c("blue", "pink", "black", "purple", "greenyellow"))

#==============================================================================
# 4. Univariate Cox and KM Analysis
#==============================================================================

# KM survival analysis for each gene
surv_results <- lapply(lactylation_genes, function(gene) {
  expr <- as.numeric(tcga_expr[gene,])
  median_expr <- median(expr, na.rm = TRUE)
  group <- ifelse(expr > median_expr, "High", "Low")
  
  surv_obj <- Surv(surv_data$OS.time, surv_data$OS)
  fit <- survfit(surv_obj ~ group)
  
  # Log-rank test
  test <- survdiff(surv_obj ~ group)
  pval <- 1 - pchisq(test$chisq, df = 1)
  
  # Cox regression
  cox_model <- coxph(surv_obj ~ expr)
  cox_coef <- summary(cox_model)$coefficients
  
  list(
    gene = gene,
    km_pval = pval,
    hr = cox_coef[2],
    hr_low = summary(cox_model)$conf.int[3],
    hr_high = summary(cox_model)$conf.int[4],
    cox_pval = cox_coef[5]
  )
})

# Filter significant genes (p < 0.05 in both KM and Cox)
sig_prognostic_genes <- sapply(surv_results, function(x) {
  x$km_pval < 0.05 && x$cox_pval < 0.05
})
prognostic_genes <- lactylation_genes[sig_prognostic_genes]

#==============================================================================
# 5. LASSO Regression (Fig S8A, S8B)
#==============================================================================

# Prepare LASSO input
x <- t(tcga_expr[prognostic_genes,])
y <- Surv(surv_data$OS.time, surv_data$OS)

# Cross-validation for lambda selection
cvfit <- cv.glmnet(x, y, family = "cox", nfolds = 10)

# Optimal lambda
best_lambda <- cvfit$lambda.min

# Fit final model
final_model <- glmnet(x, y, family = "cox", lambda = best_lambda)
coef_matrix <- as.matrix(coef(final_model))
lasso_genes <- rownames(coef_matrix)[which(coef_matrix != 0)]

#==============================================================================
# 6. LRGS Construction
#==============================================================================

# Extract coefficients
actCoef <- coef_matrix[lasso_genes,]

# Calculate risk score
riskScore <- apply(tcga_expr[lasso_genes,], 2, function(x) {
  sum(x * actCoef)
})

# Add to survival data
surv_data$riskScore <- riskScore[rownames(surv_data)]

#==============================================================================
# 7. Survival Analysis with LRGS (Fig 7C, 7D)
#==============================================================================

# Optimal cutoff
res.cut <- surv_cutpoint(surv_data, time = "OS.time", event = "OS", 
                         variables = "riskScore")
res.cat <- surv_categorize(res.cut)

# KM plot
fit <- survfit(Surv(OS.time, OS) ~ riskScore, data = res.cat)
ggsurvplot(fit, data = res.cat, pval = TRUE, risk.table = TRUE,
           title = "OS by LRGS (TCGA-SKCM)")

# Time-dependent ROC
roc_result <- timeROC(T = surv_data$OS.time,
                      delta = surv_data$OS,
                      marker = surv_data$riskScore,
                      cause = 1,
                      weighting = "marginal",
                      times = c(3, 4, 5) * 365,  # 3, 4, 5 years
                      ROC = TRUE,
                      iid = TRUE)

# Plot ROC curves
plot(roc_result, time = 3 * 365, col = "red", title = FALSE)
plot(roc_result, time = 4 * 365, col = "blue", add = TRUE)
plot(roc_result, time = 5 * 365, col = "green", add = TRUE)
legend("bottomright", legend = c("3-year", "4-year", "5-year"),
       col = c("red", "blue", "green"), lwd = 2)

#==============================================================================
# 8. Multivariate Cox Regression (Fig S8E, S8F)
#==============================================================================

# Univariate Cox
uni_cox <- coxph(Surv(OS.time, OS) ~ riskScore, data = surv_data)
summary(uni_cox)

# Multivariate Cox
multi_cox <- coxph(Surv(OS.time, OS) ~ riskScore + Age + Stage + Gender, 
                   data = surv_data)
summary(multi_cox)

# Forest plot
ggforest(multi_cox, data = surv_data)

#==============================================================================
# 9. Validation in GSE65904
#==============================================================================

# Apply LRGS to validation cohort
valid_expr <- gse65904_expr[lasso_genes,]
valid_riskScore <- apply(valid_expr, 2, function(x) sum(x * actCoef))

# KM and ROC in validation set
valid_surv <- data.frame(
  OS.time = valid_clinical$OS.time,
  OS = valid_clinical$OS,
  riskScore = valid_riskScore
)

# Validation KM
res.cut.valid <- surv_cutpoint(valid_surv, time = "OS.time", event = "OS",
                               variables = "riskScore")
res.cat.valid <- surv_categorize(res.cut.valid)
fit_valid <- survfit(Surv(OS.time, OS) ~ riskScore, data = res.cat.valid)
ggsurvplot(fit_valid, data = res.cat.valid, pval = TRUE)

saveRDS(list(
  LRGS_genes = lasso_genes,
  coefficients = actCoef,
  tcga_results = surv_data,
  validation_results = valid_surv
), "results/09_LRGS_model.rds")
