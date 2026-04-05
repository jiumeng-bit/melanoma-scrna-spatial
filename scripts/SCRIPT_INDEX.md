# Script Index

Complete guide to all analysis scripts with execution order, dependencies, and expected outputs.

## Execution Order

Scripts should be run in numerical order (01 → 11) as each script depends on outputs from previous scripts.

```
01_data_preprocessing.R 
    ↓
02_TME_heterogeneity_analysis.R
    ↓
    ├→ 03_melanocyte_subtype_analysis.R ──→ 04_lactylation_correlation_analysis.R
    │                                          ↓
    │                                       09_TCGA_LRGS_analysis.R
    ├→ 05_myeloid_cell_analysis.R ─────────→ 08_cellchat_interaction.R
    │                                          ↓
    └→ 06_T_cell_analysis.R ────────────────→ 10_ICB_therapy_analysis.R
                                                  ↓
    ┌────────────────────────────────────────→ 11_SLC25A39_validation.R
    ↓
07_spatial_transcriptomics.R (Optional)
```

## Script Details

### 01_data_preprocessing.R
**Status**: ⭐ Required first step

**Input**:
- `data/Melanoma_GSE72056.txt.gz`
- `data/Melanoma_GSE115978.rds`
- `data/Melanoma_GSE215120/*.h5`
- `data/healthy_skin_PRJCA006797.rds`
- `data/gencode_v19_gene_pos.txt`

**Output**:
- `results/01_scRNA_annotated.rds` (3.5 GB)

**Runtime**: ~2 hours
**Memory**: ~32 GB

**Key Functions**:
- Data merging and QC
- Harmony batch correction
- Cell type annotation

---

### 02_TME_heterogeneity_analysis.R
**Status**: ⭐ Required

**Input**:
- `results/01_scRNA_annotated.rds`
- `data/kegg_metabolic_pathways.rds`

**Output**:
- `results/02_scRNA_scored.rds`

**Runtime**: ~30 minutes
**Memory**: ~16 GB

**Key Functions**:
- Hypoxia scoring (AUCell)
- Lactylation scoring
- GSVA pathway analysis

---

### 03_melanocyte_subtype_analysis.R
**Status**: ⭐ Required

**Input**:
- `results/02_scRNA_scored.rds`

**Output**:
- `results/03_melanocyte_subtype.rds`
- `results/infercnv_output/` (folder)
- `results/melanocyte_counts.csv` (for cNMF)
- `results/melanocyte_scfea_input.csv` (for scFEA)

**Runtime**: ~4 hours (includes inferCNV)
**Memory**: ~32 GB

**Key Functions**:
- Melanocyte subclustering (M1-M6)
- inferCNV analysis
- scMetabolism scoring

**Note**: Requires external tools:
- cNMF (Python): `cnmf cnmf --output-dir results/cnmf/ --name melanocyte`
- scFEA (Python): `python scFEA.py --input_file results/melanocyte_scfea_input.csv`

---

### 04_lactylation_correlation_analysis.R
**Status**: ⭐ Required

**Input**:
- `results/03_melanocyte_subtype.rds`
- `results/cnmf_usage.csv` (from cNMF)

**Output**:
- `results/04_melanocyte_lactylation.rds`

**Runtime**: ~3 hours (includes WGCNA)
**Memory**: ~32 GB

**Key Functions**:
- Monocle3 pseudotime analysis
- hdWGCNA network construction
- Module-trait correlation

**Note**: Requires hdWGCNA installation:
```r
devtools::install_github('smorabit/hdWGCNA', ref='dev')
```

---

### 05_myeloid_cell_analysis.R
**Status**: ⭐ Required

**Input**:
- `results/02_scRNA_scored.rds`

**Output**:
- `results/05_myeloid_analysis.rds`
- `results/myeloid_counts_for_pyscenic.csv`

**Runtime**: ~1 hour
**Memory**: ~16 GB

**Key Functions**:
- Myeloid subclustering
- scMetabolism analysis
- Monocle3 trajectory

**Note**: Exports data for pySCENIC (Python)

---

### 06_T_cell_analysis.R
**Status**: ⭐ Required

**Input**:
- `results/02_scRNA_scored.rds`

**Output**:
- `results/06_T_cell_analysis.rds`

**Runtime**: ~1 hour
**Memory**: ~16 GB

**Key Functions**:
- T cell subclustering
- CD8_Tact characterization
- Exhaustion marker analysis

---

### 07_spatial_transcriptomics.R
**Status**: 📊 Optional (requires spatial data)

**Input**:
- `results/01_scRNA_annotated.rds` (reference)
- `data/spatial/NS1/` (10x Visium data)
- `data/spatial/SM1/` (10x Visium data)
- etc.

**Output**:
- `results/07_spatial_analysis.rds`
- `results/rctd_results/` (folder)

**Runtime**: ~2 hours per sample
**Memory**: ~16 GB

**Key Functions**:
- RCTD deconvolution
- Spatial metabolic scoring
- Distance-based analysis

---

### 08_cellchat_interaction.R
**Status**: ⭐ Required

**Input**:
- `results/02_scRNA_scored.rds`
- `results/03_melanocyte_subtype.rds`
- `results/05_myeloid_analysis.rds`
- `results/06_T_cell_analysis.rds`

**Output**:
- `results/08_cellchat_analysis.rds`

**Runtime**: ~2 hours
**Memory**: ~32 GB

**Key Functions**:
- CellChat network inference
- MIF pathway analysis
- CXCL16-CXCR6 axis
- CCL5-CCR5 axis

---

### 09_TCGA_LRGS_analysis.R
**Status**: ⭐ Required

**Input**:
- TCGA/GTEx data (downloaded via TCGAbiolinks)
- `results/04_melanocyte_lactylation.rds` (WGCNA modules)

**Output**:
- `results/09_LRGS_model.rds`

**Runtime**: ~1 hour
**Memory**: ~16 GB

**Key Functions**:
- Differential expression
- LASSO regression
- Survival analysis
- ROC curves

---

### 10_ICB_therapy_analysis.R
**Status**: 📊 Optional (requires ICB cohorts)

**Input**:
- `data/GSE120575_seurat.rds`
- IMvigor210 package data
- `results/09_LRGS_model.rds`

**Output**:
- `results/10_ICB_analysis.rds`

**Runtime**: ~30 minutes
**Memory**: ~16 GB

**Key Functions**:
- Pre/On-treatment comparison
- Responder vs non-responder analysis
- LRGS prediction

---

### 11_SLC25A39_validation.R
**Status**: 📊 Optional (requires experimental data)

**Input**:
- `data/experiments/SLC25A39_qPCR_results.csv`
- `data/experiments/SLC25A39_CCK8_results.csv`
- `data/experiments/SLC25A39_transwell_results.csv`
- `data/experiments/SLC25A39_wound_healing_results.csv`
- `data/experiments/SLC25A39_WB_results.csv`

**Output**:
- `results/11_SLC25A39_validation.rds`
- Publication figures

**Runtime**: ~15 minutes
**Memory**: ~4 GB

**Key Functions**:
- qPCR analysis
- Proliferation curves
- Migration/invasion quantification
- Western blot densitometry

---

## Parallel Execution

Some scripts can be run in parallel after completing dependencies:

```
After 02_TME_heterogeneity_analysis.R:
├── 03_melanocyte_subtype_analysis.R ──┬── 04_lactylation_correlation_analysis.R
│                                      └── 08_cellchat_interaction.R
├── 05_myeloid_cell_analysis.R ────────┘
└── 06_T_cell_analysis.R ──────────────┘
```

## Error Recovery

If a script fails:

1. Check the error message
2. Verify input files exist
3. Check memory availability
4. Re-run from the failed script (previous outputs are preserved)

## Resource Monitoring

Monitor resource usage during execution:

```bash
# In a separate terminal
watch -n 5 'free -h && df -h'
```

## Progress Tracking

Create a progress file:

```bash
# Start analysis
echo "Starting analysis: $(date)" > analysis_progress.log

# After each script
echo "Completed 01: $(date)" >> analysis_progress.log
```
