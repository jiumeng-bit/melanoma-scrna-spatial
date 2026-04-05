#==============================================================================
# Script 7: Spatial Transcriptomics Analysis
# Description: Spatial mapping, deconvolution, metabolic landscape
# Corresponds to: Results 3.6
#==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(SPOTlight)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(dplyr)
})

#==============================================================================
# 1. Load Spatial Data (Fig 6A)
#==============================================================================

# GSE225475 - Healthy skin
# GSE179572, GSE203612 - Melanoma

# Example for one sample
spatial_data <- Load10X_Spatial(
  data.dir = "data/spatial/NS1",
  filename = "filtered_feature_bc_matrix.h5"
)

# QC
spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")
spatial_data <- subset(spatial_data, 
                       subset = nFeature_Spatial > 5 & 
                         nFeature_Spatial < 300 & 
                         percent.mt < 20)

# SCTransform normalization
spatial_data <- SCTransform(spatial_data, assay = "Spatial", verbose = FALSE)

#==============================================================================
# 2. RCTD Deconvolution (Fig 6A, 6B)
#==============================================================================

# Prepare reference from scRNA-seq
scRNA <- readRDS("results/01_scRNA_annotated.rds")

# Get marker genes
markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

# Create RCTD reference
ref_counts <- GetAssayData(scRNA, layer = "counts")
ref_cell_types <- scRNA$celltype
names(ref_cell_types) <- colnames(scRNA)

reference <- Reference(ref_counts, ref_cell_types)

# Create SpatialRNA object
spatial_counts <- GetAssayData(spatial_data, layer = "counts")
spatial_coords <- GetTissueCoordinates(spatial_data)
spatial_coords <- spatial_coords[, c("imagerow", "imagecol")]

coords <- data.frame(
  x = spatial_coords$imagecol,
  y = spatial_coords$imagerow
)
rownames(coords) <- colnames(spatial_counts)

puck <- SpatialRNA(coords, spatial_counts)

# Run RCTD
my_rctd <- create.RCTD(puck, reference, max_cores = 8)
my_rctd <- run.RCTD(my_rctd, doublet_mode = "full")

# Get results
results <- my_rctd@results
norm_weights <- normalize_weights(results$weights)
cell_type_names <- my_rctd@cell_type_info$info[[2]]
spatial_data@meta.data <- cbind(spatial_data@meta.data, norm_weights)

#==============================================================================
# 3. Spatial Visualization (Fig 6C)
#==============================================================================

# Plot cell type proportions
p1 <- SpatialFeaturePlot(spatial_data, features = "Melanocyte", 
                         pt.size.factor = 1.5, alpha = c(0.1, 1))
p2 <- SpatialFeaturePlot(spatial_data, features = "CD8_Tact",
                         pt.size.factor = 1.5, alpha = c(0.1, 1))

#==============================================================================
# 4. Spatial Metabolic Analysis (Fig 6D, 6E, 6F, 6G)
#==============================================================================

# Calculate lactylation score in spatial data
lactylation_genes <- c("HIF1A", "LDHA", "ENO1", "GAPDH", "PGK1", "PGAM1", "PKM", "HK2")
spatial_data <- AddModuleScore(spatial_data, features = list(lactylation_genes),
                               name = "lactylation_score")

# Plot lactylation
p3 <- SpatialFeaturePlot(spatial_data, features = "lactylation_score1",
                         pt.size.factor = 1.5)

# Distance from tumor center analysis
# Identify tumor center based on melanocyte density
tumor_spots <- which(spatial_data$Melanocyte > 0.3)
coords <- GetTissueCoordinates(spatial_data)
tumor_center <- colMeans(coords[tumor_spots, c("imagerow", "imagecol")])

# Calculate distance from center
spatial_data$dist_from_center <- sqrt(
  (coords$imagerow - tumor_center["imagerow"])^2 +
    (coords$imagecol - tumor_center["imagecol"])^2
)

# Correlation between lactylation and distance (Fig 6E)
distance_lactylation <- data.frame(
  distance = spatial_data$dist_from_center,
  lactylation = spatial_data$lactylation_score1
)

p4 <- ggplot(distance_lactylation, aes(x = distance, y = lactylation)) +
  geom_smooth(method = "loess", color = "red") +
  geom_point(alpha = 0.3) +
  theme_bw() +
  labs(title = "Lactylation vs Distance from Tumor Center")

saveRDS(spatial_data, "results/07_spatial_analysis.rds")
