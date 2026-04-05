# Single-Cell and Spatial Transcriptomic Analysis Reveals Lactylation-Driven Metabolic Reprogramming in Skin Cutaneous Melanoma


This repository contains the complete analysis pipeline for the manuscript investigating lactylation-driven metabolic reprogramming in skin cutaneous melanoma (SKCM) through integration of single-cell RNA sequencing, spatial transcriptomics, and bulk RNA-seq data.

## 📋 Overview

This study reveals how lactylation and glycolytic metabolic reprogramming drive melanoma progression through:

- **Multi-omics integration**: 4 scRNA-seq datasets + 3 spatial transcriptomics datasets + TCGA/GTEx bulk RNA-seq
- **Key discoveries**:
  - **M5_ASTN2**: Highly malignant melanocyte subtype with peak glycolytic activity
  - **Mφ_CCL5**: Glycolytic macrophage recruiting CD8+ T cells via CXCL16-CXCR6 axis
  - **LRGS**: Novel 8-gene prognostic signature for survival and immunotherapy response
  - **SLC25A39**: Functionally validated as key lactylation regulator

## 🗂️ Repository Structure

```
├── scripts/
│   ├── 01_data_preprocessing.R              # Data integration & QC
│   ├── 02_TME_heterogeneity_analysis.R      # TME characterization
│   ├── 03_melanocyte_subtype_analysis.R     # Melanocyte subtypes (M1-M6)
│   ├── 04_lactylation_correlation_analysis.R # Trajectory & WGCNA
│   ├── 05_myeloid_cell_analysis.R           # Myeloid reprogramming
│   ├── 06_T_cell_analysis.R                 # CD8+ T cell exhaustion
│   ├── 07_spatial_transcriptomics.R         # Spatial metabolic landscape
│   ├── 08_cellchat_interaction.R            # Cell-cell communication
│   ├── 09_TCGA_LRGS_analysis.R              # Prognostic model construction
│   ├── 10_ICB_therapy_analysis.R            # Immunotherapy response
```

## 🚀 Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/[username]/melanoma-scrna-spatial.git
cd melanoma-scrna-spatial
```

### 2. Install Dependencies

```r
# Install CRAN packages
install.packages(c("Seurat", "harmony", "survival", "survminer", "glmnet"))

# Install Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "infercnv", "DESeq2", "sva"))

# Install GitHub packages
remotes::install_github("genecell/COSG")
remotes::install_github("sqjin/CellChat")
```



## 📊 Key Results

### 1. Melanocyte Heterogeneity
- **6 subtypes identified**: M1_SPARCL1, M2_CRABP1, M3_ABCB5, M4_ISG15, **M5_ASTN2**, M6_MKI67
- **M5 characteristics**: Highest CNV score, glycolytic metaprogram (LDHA+), worst prognosis

### 2. Immune Metabolic Reprogramming
- **Mφ_CCL5**: Highly glycolytic, immunosuppressive (ARG1+, CD274+)
- **CD8_Tact**: Transitional exhaustion state, CXCR6+ for recruitment

### 3. Spatial Metabolic Landscape
- **Tumor core**: Highest lactylation, inhibits immune infiltration
- **Tumor edge**: Mφ_CCL5 and CD8_Tact co-localization

### 4. LRGS Prognostic Model
- **8 genes**: TIMM50, SLC25A39, CMTM8, ICAM1, COMMD3, SPATA13, EGR3, HLA-DMA
- **Performance**: AUC 0.69-0.71 for 3/4/5-year survival
- **Validation**: TCGA + GSE65904

### 5. Immunotherapy Prediction
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

## 📝 Documentation

- **[ANALYSIS_OVERVIEW.md](docs/ANALYSIS_OVERVIEW.md)**: Detailed method descriptions
- **[QUICKSTART.md](docs/QUICKSTART.md)**: Step-by-step setup guide
- **[SCRIPT_INDEX.md](scripts/SCRIPT_INDEX.md)**: Script execution order and dependencies

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## ⚠️ Notes

- All personal file paths have been anonymized for GitHub publication
- Please update paths in scripts according to your system configuration
- Some analyses require Python tools (cNMF, scFEA, pySCENIC)


