# Single-Cell and Spatial Transcriptomic Analysis Reveals Lactylation-Driven Metabolic Reprogramming in Skin Cutaneous Melanoma

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

This repository contains the complete bioinformatics analysis pipeline for investigating lactylation-driven metabolic reprogramming in skin cutaneous melanoma (SKCM) through integration of single-cell RNA sequencing, spatial transcriptomics, and bulk RNA-seq data.

**Key discoveries**:
- **M5_ASTN2**: Highly malignant melanocyte subtype with peak glycolytic activity
- **Mφ_CCL5**: Glycolytic macrophage recruiting CD8+ T cells via CXCL16-CXCR6 axis
- **LRGS**: Novel 8-gene prognostic signature for survival and immunotherapy response

## 🗂️ Repository Structure

```
.
├── scripts/                    # 10 R analysis scripts
│   ├── 01_data_preprocessing.R              # Data integration & QC
│   ├── 02_TME_heterogeneity_analysis.R      # TME characterization
│   ├── 03_melanocyte_subtype_analysis.R     # Melanocyte subtypes (M1-M6)
│   ├── 04_lactylation_correlation_analysis.R # Trajectory & WGCNA
│   ├── 05_myeloid_cell_analysis.R           # Myeloid reprogramming
│   ├── 06_T_cell_analysis.R                 # CD8+ T cell exhaustion
│   ├── 07_spatial_transcriptomics.R         # Spatial metabolic landscape
│   ├── 08_cellchat_interaction.R            # Cell-cell communication
│   ├── 09_TCGA_LRGS_analysis.R              # Prognostic model construction
│   └── 10_ICB_therapy_analysis.R            # Immunotherapy response
│
├── data/                       # Data directory (user-provided)
├── results/                    # Output directory
└── README.md                   # This file
```

## 📊 Key Results

### 1. Melanocyte Heterogeneity (Fig 2)
- **6 subtypes identified**: M1_SPARCL1, M2_CRABP1, M3_ABCB5, M4_ISG15, **M5_ASTN2**, M6_MKI67
- **M5 characteristics**: Highest CNV score, glycolytic metaprogram (LDHA+), worst prognosis

### 2. Immune Metabolic Reprogramming (Fig 4-5)
- **Mφ_CCL5**: Highly glycolytic, immunosuppressive (ARG1+, CD274+)
- **CD8_Tact**: Transitional exhaustion state, CXCR6+ for recruitment

### 3. Spatial Metabolic Landscape (Fig 6)
- **Tumor core**: Highest lactylation, inhibits immune infiltration
- **Tumor edge**: Mφ_CCL5 and CD8_Tact co-localization

### 4. LRGS Prognostic Model (Fig 7)
- **8 genes**: TIMM50, SLC25A39, CMTM8, ICAM1, COMMD3, SPATA13, EGR3, HLA-DMA
- **Performance**: AUC 0.69-0.71 for 3/4/5-year survival
- **Validation**: TCGA + GSE65904

### 5. Immunotherapy Prediction (Fig 9)
- High lactylation correlates with poor ICB response
- SLC25A39 decreases post-treatment in responders

## 📚 Data Availability

| Data Type | Accession | Sample Size |
|-----------|-----------|-------------|
| scRNA-seq | GSE72056 | 19 samples |
| scRNA-seq | GSE115978 | 16 samples |
| scRNA-seq | GSE215120 | 4 samples |
| scRNA-seq | PRJCA006797 | 5 healthy |
| Spatial | GSE225475 | 2 healthy |
| Spatial | GSE179572 | 3 primary |
| Spatial | GSE203612 | 1 metastasis |
| Bulk RNA | TCGA-SKCM | 470 samples |
| Bulk RNA | GTEx | Normal skin |
| Bulk RNA | GSE65904 | 210 samples |
| ICB scRNA | GSE120575 | Pre/On treatment |
| ICB Bulk | IMvigor210 | Urothelial cancer |

## 🔬 Methods Highlights

### Single-Cell Analysis
- **Batch correction**: Harmony (sample-level)
- **Clustering**: Seurat (resolution 0.1-0.8)
- **Annotation**: COSG + manual curation

### Advanced Analyses
- **CNV**: inferCNV with M2 as reference
- **Metaprograms**: cNMF (7 programs)
- **Metabolism**: scMetabolism + scFEA
- **Trajectory**: Monocle3 + hdWGCNA

### Spatial Analysis
- **Deconvolution**: RCTD (spacexr)
- **Cell mapping**: 9 cell types annotated
- **Gradient analysis**: Distance from tumor center

### Statistical Methods
- **Survival**: Kaplan-Meier, Cox regression
- **Feature selection**: LASSO regression
- **ROC analysis**: timeROC for time-dependent AUC
