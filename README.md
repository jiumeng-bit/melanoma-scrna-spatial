# Single-Cell and Spatial Transcriptomic Analysis Reveals Lactylation-Driven Metabolic Reprogramming in Skin Cutaneous Melanoma

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete bioinformatics analysis pipeline for the manuscript investigating lactylation-driven metabolic reprogramming in skin cutaneous melanoma (SKCM) through integration of single-cell RNA sequencing, spatial transcriptomics, and bulk RNA-seq data.

## 📋 Overview

This study reveals how lactylation and glycolytic metabolic reprogramming drive melanoma progression through:

- **Multi-omics integration**: 4 scRNA-seq datasets + 3 spatial transcriptomics datasets + TCGA/GTEx bulk RNA-seq
- **56,168 single cells** across 44 samples (39 melanoma, 5 healthy skin)
- **Key discoveries**:
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
├── docs/                       # Documentation
│   ├── ANALYSIS_OVERVIEW.md     # Detailed workflow documentation
│   └── QUICKSTART.md           # Getting started guide
│
├── data/                       # Data directory (user-provided)
├── results/                    # Output directory
├── requirements.txt            # R package dependencies
└── README.md                   # This file
```

## 🚀 Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/jiumeng-bit/melanoma-scrna-spatial.git
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

See `requirements.txt` for complete list.

### 3. Run Analysis

```r
# Run complete pipeline
for(i in 1:10) {
  source(sprintf("scripts/%02d_*.R", i))
}

# Or run specific analysis
source("scripts/03_melanocyte_subtype_analysis.R")
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

## 📝 Documentation

- **[ANALYSIS_OVERVIEW.md](docs/ANALYSIS_OVERVIEW.md)**: Detailed method descriptions
- **[QUICKSTART.md](docs/QUICKSTART.md)**: Step-by-step setup guide
- **[SCRIPT_INDEX.md](scripts/SCRIPT_INDEX.md)**: Complete script execution guide

## 🖥️ System Requirements

- **OS**: Linux/macOS/Windows
- **R**: >= 4.0.0
- **RAM**: 32 GB recommended (64 GB for full pipeline)
- **CPU**: Multi-core recommended
- **Disk**: 50 GB free space

## ⚡ Runtime

| Script | Runtime | Memory |
|--------|---------|--------|
| 01 | 2 hours | 32 GB |
| 02 | 30 min | 16 GB |
| 03 | 4 hours | 32 GB |
| 04-06 | 1-3 hours each | 16-32 GB |
| 07 | 2 hours/sample | 16 GB |
| 08 | 2 hours | 32 GB |
| 09-10 | 15-60 min | 4-16 GB |

**Total**: ~16-20 hours for complete pipeline

## 📖 Citation

```bibtex
@article{yu2024lactylation,
  title={Single-Cell and Spatial Transcriptomic Analysis Reveals Lactylation-Driven Metabolic Reprogramming in Skin Cutaneous Melanoma},
  author={Yu, Yongkai and Wang, Wenjing and Wang, Yidan and Pei, Tongxin and Lu, Jiawei and Feng, Yifei and Cao, Xuechen and Lu, Yan},
  journal={Translational Research},
  year={2024},
  publisher={Elsevier}
}
```

## 👥 Contributors

- **Yongkai Yu**: Conceptualization, Data curation, Formal analysis, Visualization
- **Wenjing Wang**: Validation, Methodology
- **Yidan Wang**: Formal analysis, Visualization
- **Yan Lu**: Supervision, Funding acquisition, Conceptualization

## 📧 Contact

**Corresponding Author**: Yan Lu  
**Email**: luyan6289@163.com  
**Institution**: Department of Dermatology, The First Affiliated Hospital of Nanjing Medical University

## 🙏 Acknowledgments

- National Natural Science Foundation of China (No. 82273549)
- Natural Science Foundation of Jiangsu Province (SBK2022023090)
- Jiangsu Provincial Medical Key Discipline Cultivation Unit (JSDW202228)

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## ⚠️ Notes

- All personal file paths have been anonymized for GitHub publication
- Please update paths in scripts according to your system configuration
- Some analyses require Python tools (cNMF, scFEA, pySCENIC)

---

**Last Updated**: April 2024  
**Version**: 1.0.0
