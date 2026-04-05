# Detailed Analysis Overview

## Manuscript Structure Correspondence

| Script | Manuscript Section | Key Figures | Description |
|--------|-------------------|-------------|-------------|
| 01_data_preprocessing.R | Methods 2.1 | Fig 1A, S1 | Data integration and QC |
| 02_TME_heterogeneity_analysis.R | Results 3.1 | Fig 1B-F, S1 | TME characterization |
| 03_melanocyte_subtype_analysis.R | Results 3.2 | Fig 2, S2 | Melanocyte heterogeneity |
| 04_lactylation_correlation_analysis.R | Results 3.3 | Fig 3, S3 | Lactylation dynamics |
| 05_myeloid_cell_analysis.R | Results 3.4 | Fig 4, S4 | Myeloid reprogramming |
| 06_T_cell_analysis.R | Results 3.5 | Fig 5, S5 | T cell exhaustion |
| 07_spatial_transcriptomics.R | Results 3.6 | Fig 6, S6, S7 | Spatial metabolism |
| 08_cellchat_interaction.R | Results 3.4-3.5 | Fig 4F-G, 5I-J, S4C-D | Cell communication |
| 09_TCGA_LRGS_analysis.R | Results 3.7 | Fig 7, S8 | Prognostic model |
| 10_ICB_therapy_analysis.R | Results 3.8 | Fig 9, S9 | Immunotherapy |

## Detailed Workflow

### 1. Data Preprocessing (01_data_preprocessing.R)

**Input**: Raw count matrices from 4 datasets
- GSE72056: 19 samples (7 primary, 12 metastatic)
- GSE115978: 16 samples (4 primary, 12 metastatic)
- GSE215120: 4 samples (primary melanoma)
- PRJCA006797: 5 healthy skin samples

**Output**: Integrated Seurat object with annotations

**Key Steps**:
1. Quality control (nCount_RNA: 500-10000, nFeature_RNA: 200-5000, MT% < 20)
2. Harmony batch correction by sample
3. Clustering (resolution = 0.8)
4. Cell type annotation using COSG + manual curation

**Cell Types Identified**:
- Melanocyte (PMEL, TYRP1, MLANA)
- Keratinocyte (KRT1, KRT5, KRT10)
- Fibroblast (COL1A1, COL1A2, DCN)
- Endothelial (PECAM1, VWF)
- T cell (CD3D, CD3E)
- B cell (CD79A, MS4A1)
- Myeloid (LYZ, CD14, CD68)
- Langerhans (CD1A, CD207)

### 2. TME Heterogeneity (02_TME_heterogeneity_analysis.R)

**Key Analyses**:
1. **Cell Proportion Analysis**: Bar charts showing changes between healthy and tumor
2. **K-means Clustering**: Identified 4 microenvironment subtypes (Endothelialhi, T cellhi, Myeloidhi, B cellhi)
3. **Hypoxia Scoring**: Using HIF1A, EPAS1, VEGFA, etc.
4. **Lactylation Scoring**: Using 8 lactylation-related genes
5. **Glycolysis Scoring**: GSVA with KEGG pathway

**Key Findings**:
- Melanocytes show highest hypoxia and lactylation in tumor
- HIF1A+ melanocytes have enhanced lactylation
- Glycolysis upregulated specifically in tumor melanocytes

### 3. Melanocyte Subtypes (03_melanocyte_subtype_analysis.R)

**Subtypes Identified** (resolution = 0.1):
- M1_SPARCL1: Low CNV, normal-like
- M2_CRABP1: Healthy-enriched, reference for inferCNV
- M3_ABCB5: Moderate CNV
- M4_ISG15: Immune-related
- M5_ASTN2: **Highest CNV, glycolytic, most malignant**
- M6_MKI67: Proliferating

**Analyses**:
1. **inferCNV**: CNV quantification with M2 as reference
2. **cNMF**: 7 metaprograms identified, C1 = Glycolytic (LDHA-enriched)
3. **scMetabolism**: KEGG pathway activities
4. **scFEA**: Metabolic flux analysis (pyruvate to lactate shift in M5)

**Key Findings**:
- M5 has highest CNV score and glycolytic activity
- Glycolytic metaprogram correlates with poor prognosis
- ALDOA identified as hub gene for glycolytic flux

### 4. Lactylation Trajectory (04_lactylation_correlation_analysis.R)

**Analyses**:
1. **Pathway Correlation**: Lactylation vs Glycolysis, mTOR, TGF-β, EMT
2. **Monocle3 Pseudotime**: M2 → M1/M3-M6 trajectory, M5 at terminal stage
3. **hdWGCNA**: Scale-free network (soft threshold = 5), 12 modules identified
   - Lactylation-related: blue, pink, black, purple, greenyellow (cor > 0.4)

**Key Findings**:
- Strong positive correlation between lactylation and glycolysis
- Lactylation increases along pseudotime
- ENO1 and ALDH1A1 enriched in late pseudotime

### 5. Myeloid Analysis (05_myeloid_cell_analysis.R)

**Subtypes**:
- cDC1_CLEC9A: Conventional DC1
- cDC2_CD1C: Conventional DC2
- Mono_MKI67: Cycling monocytes
- Mφ_C1QB: Resident macrophages
- Mφ_CXCL8: Inflammatory macrophages
- **Mφ_CCL5: Glycolytic, immunoregulatory**

**Key Findings**:
- Mφ_CCL5 shows most active glycolysis
- Expresses ARG1, MRC1, CD274, CX3CR1 (immunosuppressive markers)
- Trajectory shows increasing glycolysis/OXPHOS along differentiation

### 6. T Cell Analysis (06_T_cell_analysis.R)

**Subtypes**:
- CD4_Treg: Regulatory T cells
- CD8_Tcyt: Cytotoxic T cells
- **CD8_Tact: Activated, glycolytic, transitional state**
- CD8_Tpro: Proliferating
- CD8_Tifn: IFN-producing
- CD8_Tex: Exhausted

**Key Findings**:
- CD8_Tact enriched in tumors, shows dual activation/exhaustion phenotype
- High glycolysis correlates with PDCD1, TIGIT, HAVCR2 expression
- Expresses CCR5, CXCR6 (recruitment markers)

### 7. Spatial Transcriptomics (07_spatial_transcriptomics.R)

**Samples**:
- 2 healthy skin (GSE225475)
- 3 primary SKCM (10x Genomics)
- 1 brain metastasis (GSE203612)

**Analyses**:
1. **RCTD Deconvolution**: Cell type proportions per spot
2. **Spatial Mapping**: M5 enrichment in tumor core
3. **Mφ_CCL5-CD8_Tact Co-localization**: At tumor edge
4. **Distance Analysis**: Lactylation decreases from tumor center

**Key Findings**:
- M5 concentrated in tumor core
- Mφ_CCL5 and CD8_Tact co-localize at tumor edge
- High lactylation environment inhibits immune infiltration

### 8. Cell-Cell Communication (08_cellchat_interaction.R)

**Interactions Analyzed**:
1. **MIF Pathway**: M5 → Mφ_CCL5 (MIF-CD74+CXCR4, MIF-CD74+CD44)
2. **CXCL Pathway**: Mφ_CCL5 → CD8_Tact (CXCL16-CXCR6)
3. **CCL Pathway**: Mφ_CCL5 → CD8_Tact (CCL5-CCR5)

**Key Findings**:
- M5-Mφ_CCL5 interaction strongest in tumor
- MIF pathway enriched in tumor vs healthy
- Mφ_CCL5 recruits CD8_Tact via CXCL16-CXCR6 axis

### 9. TCGA Analysis and LRGS (09_TCGA_LRGS_analysis.R)

**Cohorts**:
- TCGA-SKCM (training): 470 samples
- GSE65904 (validation): 210 samples

**Pipeline**:
1. Differential expression (DESeq2)
2. Univariate Cox/KM survival analysis
3. **LASSO regression**: Selected 8 genes (TIMM50, SLC25A39, CMTM8, ICAM1, COMMD3, SPATA13, EGR3, HLA-DMA)
4. LRGS construction: Weighted sum of 8 genes
5. Time-dependent ROC (3/4/5-year)
6. Multivariate Cox regression

**Key Findings**:
- LRGS predicts OS (AUC: 0.69-0.71)
- Independent prognostic factor (HR > 1, p < 0.001)
- SLC25A39 only gene remaining significant after multivariate correction

### 10. ICB Therapy Analysis (10_ICB_therapy_analysis.R)

**Cohorts**:
- GSE120575 (scRNA-seq): Pre/On treatment, responders/non-responders
- IMvigor210 (bulk RNA): Atezolizumab-treated urothelial cancer (validation)

**Key Findings**:
- Non-responders have higher baseline lactylation
- Responders show increased CD8_Tact post-treatment
- SLC25A39 decreases after treatment in responders
- High LRGS correlates with poor ICB response (84% PD/SD)

## Dependencies

### R Packages

```r
# Core
library(Seurat)           # Single-cell analysis
library(harmony)          # Batch correction
library(SingleCellExperiment)

# Analysis
library(COSG)             # Marker identification
library(infercnv)         # CNV analysis
library(scMetabolism)     # Metabolic pathway scoring
library(Monocle3)         # Trajectory analysis
library(hdWGCNA)          # Co-expression network
library(CellChat)         # Cell-cell communication
library(spacexr)          # Spatial deconvolution (RCTD)

# Survival and statistics
library(survival)
library(survminer)
library(glmnet)           # LASSO
library(timeROC)
library(DESeq2)

# Visualization
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
```

### Python Packages (for cNMF, scFEA, pySCENIC)

```bash
# cNMF
pip install cnmf

# scFEA
# Follow installation at: https://github.com/changwn/scFEA

# pySCENIC
pip install pyscenic
```

## Computational Requirements

- **Memory**: 64GB RAM recommended for large-scale integration
- **CPU**: Multi-core processing (8+ cores)
- **Runtime**: Full pipeline ~16-20 hours depending on dataset size

## Notes

1. **Random Seeds**: Set for reproducibility in clustering and LASSO
2. **Batch Correction**: Harmony used consistently across analyses
3. **Statistical Tests**: Wilcoxon for pairwise, ANOVA/Kruskal-Wallis for multiple groups
4. **Multiple Testing**: Benjamini-Hochberg FDR correction applied
