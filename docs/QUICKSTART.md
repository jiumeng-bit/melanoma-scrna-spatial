# Quick Start Guide

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/jiumeng-bit/melanoma-scrna-spatial.git
cd melanoma-scrna-spatial
```

### 2. Install R Dependencies

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
for(i in 1:11) {
  source(sprintf("scripts/%02d_*.R", i))
}

# Or run specific analysis
source("scripts/03_melanocyte_subtype_analysis.R")
```

## Data Preparation

### Required Data Structure

```
data/
├── raw/
│   ├── Melanoma_GSE72056.txt.gz
│   ├── Melanoma_GSE115978.rds
│   ├── Melanoma_GSE215120/*.h5
│   ├── healthy_skin_PRJCA006797.rds
│   └── gencode_v19_gene_pos.txt
├── processed/          # Output directory
├── spatial/            # Spatial data (optional)
└── experiments/        # Validation data (optional)
```

### Download Public Data

```r
# Use GEOquery for GEO datasets
library(GEOquery)
gse72056 <- getGEO("GSE72056")
gse115978 <- getGEO("GSE115978")

# TCGA data via TCGAbiolinks
library(TCGAbiolinks)
# See script 09 for download code
```

## Running the Pipeline

### Option 1: Sequential Execution (Recommended for First Run)

```r
# Set working directory
setwd("path/to/melanoma-scrna-spatial")

# Run scripts in order
source("scripts/01_data_preprocessing.R")
source("scripts/02_TME_heterogeneity_analysis.R")
source("scripts/03_melanocyte_subtype_analysis.R")
source("scripts/04_lactylation_correlation_analysis.R")
source("scripts/05_myeloid_cell_analysis.R")
source("scripts/06_T_cell_analysis.R")
source("scripts/07_spatial_transcriptomics.R")      # Optional
source("scripts/08_cellchat_interaction.R")
source("scripts/09_TCGA_LRGS_analysis.R")
source("scripts/10_ICB_therapy_analysis.R")         # Optional
source("scripts/11_SLC25A39_validation.R")          # Optional
```

### Option 2: Parallel Execution (After Script 02)

Scripts 03, 05, and 06 can run in parallel after completing script 02:

```bash
# Terminal 1
Rscript scripts/03_melanocyte_subtype_analysis.R

# Terminal 2
Rscript scripts/05_myeloid_cell_analysis.R

# Terminal 3
Rscript scripts/06_T_cell_analysis.R
```

### Option 3: Run Specific Analysis

```r
# Load pre-computed data
scRNA <- readRDS("results/02_scRNA_scored.rds")

# Run only melanocyte analysis
source("scripts/03_melanocyte_subtype_analysis.R")
```

## Expected Runtime

| Script | Runtime | Memory | Output Size |
|--------|---------|--------|-------------|
| 01 | 2 hours | 32 GB | 3.5 GB |
| 02 | 30 min | 16 GB | 3.5 GB |
| 03 | 4 hours | 32 GB | 2 GB |
| 04 | 3 hours | 32 GB | 1.5 GB |
| 05 | 1 hour | 16 GB | 1 GB |
| 06 | 1 hour | 16 GB | 1 GB |
| 07 | 2 hours/sample | 16 GB | 1 GB/sample |
| 08 | 2 hours | 32 GB | 500 MB |
| 09 | 1 hour | 16 GB | 200 MB |
| 10 | 30 min | 16 GB | 100 MB |
| 11 | 15 min | 4 GB | 50 MB |

**Total**: ~18-24 hours for complete pipeline

## Troubleshooting

### Memory Issues

```r
# Increase memory limit (if supported)
options(future.globals.maxSize = 40000 * 1024^2)  # 40 GB

# Use smaller number of workers
plan(multisession, workers = 4)
```

### Package Installation Failures

```r
# If CellChat installation fails
remotes::install_github("sqjin/CellChat", force = TRUE)

# If spacexr fails (requires C++ compiler)
# Install Rtools (Windows) or Xcode (Mac) first
```

### Missing Data Files

If you don't have access to all datasets, you can:
1. Skip scripts requiring missing data
2. Use the example data from SeuratData:
```r
SeuratData::InstallData("ifnb")
ifnb <- SeuratData::LoadData("ifnb")
```

## Output Files

Results are saved in `results/` directory:

```
results/
├── 01_scRNA_annotated.rds
├── 02_scRNA_scored.rds
├── 03_melanocyte_subtype.rds
├── 04_melanocyte_lactylation.rds
├── 05_myeloid_analysis.rds
├── 06_T_cell_analysis.rds
├── 07_spatial_analysis.rds
├── 08_cellchat_analysis.rds
├── 09_LRGS_model.rds
├── 10_ICB_analysis.rds
├── 11_SLC25A39_validation.rds
├── figures/
│   ├── figure1_TME_overview.pdf
│   ├── figure2_melanocyte_subtypes.pdf
│   └── ...
└── tables/
    ├── table1_sample_info.csv
    └── ...
```

## Getting Help

1. Check `docs/ANALYSIS_OVERVIEW.md` for detailed method descriptions
2. Check `scripts/SCRIPT_INDEX.md` for script-specific information
3. Open an issue on GitHub with error messages

## Next Steps

After running the pipeline:
1. Explore results in `results/figures/`
2. Generate publication figures
3. Reproduce manuscript tables
