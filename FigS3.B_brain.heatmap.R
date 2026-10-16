# =============================================================================
# DAR Accessibility Heatmap — Brain
# =============================================================================
# Description: Generates a pseudobulk chromatin accessibility heatmap for
#              differentially accessible regions (DARs) across brain cell
#              types at 3 dpf and adult stages in zebrafish brain.
#
# Input:
#   - brain.integrated.rds
#   - "brain.merged_df.csv"  : Output file from Fig2.brain.DARs
#   - "heatmap.USE.txt"      : Cluster color palette (one hex/name
#                              per line, ordered to match clusters)
#
# Output:
#   - p5 : pheatmap object (row-scaled pseudobulk accessibility)
#
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(Seurat)
library(Signac)
library(dplyr)
library(Matrix)
library(pheatmap)


# -----------------------------------------------------------------------------
# Load object
# -----------------------------------------------------------------------------

brain <- readRDS("brain.integrated.rds")


# -----------------------------------------------------------------------------
# 1. Define exclusions
# -----------------------------------------------------------------------------

exclude_clusters <- c(
  "Basal Forebrain_C4",
  "Diencephalon_C5",
  "Cranial Ganglia_C25",
  "Cranial Motor Neurons_C26",
  "Granule Cells_C0",
  "Granule Cells_C18",
  "Granule Cells_C23",
  "Granule Cells_C36",
  "Granule Cells_C6",
  "Hypothalamic Progenitors_C38",
  "Neurons Glut+_C40",           
  "Ventral Habenula_C46"
)

# -----------------------------------------------------------------------------
# 2. Subset to 3 dpf and adult cells
# -----------------------------------------------------------------------------

brain_plot <- subset(
  brain,
  subset =
    dataset %in% c("3 dpf", "adult") &
    !(cluster.anno %in% exclude_clusters)   # note: brain uses cluster.anno
)

# Create a combined cluster_stage label (e.g. "RGCs_C4_adult")
brain_plot$cluster_stage <- paste(
  brain_plot$cluster.anno,
  brain_plot$dataset,
  sep = "_"
)


# -----------------------------------------------------------------------------
# 3. Select DAR peaks (exclude 21 dpf stage and excluded clusters)
# -----------------------------------------------------------------------------

# Load DAR results
merged_df <- read.csv("brain.merged_df.csv")

# Deduplicate: keep one row per cluster-peak combination (highest log2FC)
final_peaks <- merged_df %>%
  group_by(cluster, peak) %>%
  slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stage = sub("_vs_.*", "", comparison)
  )

peaks_filtered <- final_peaks %>%
  filter(
    stage != "21 dpf",
    !cluster %in% exclude_clusters
  )

dar_peaks <- unique(peaks_filtered$peak)


# -----------------------------------------------------------------------------
# 4. Extract accessibility matrix for DAR peaks
# -----------------------------------------------------------------------------

DefaultAssay(brain_plot) <- "ATAC"

peak_mat <- GetAssayData(
  brain_plot,
  assay = "ATAC",
  layer = "data"
)[dar_peaks, ]


# -----------------------------------------------------------------------------
# 5. Compute pseudobulk accessibility (mean per cluster_stage group)
# -----------------------------------------------------------------------------

groups       <- brain_plot$cluster_stage
group_levels <- unique(groups)

avg_mat <- do.call(
  cbind,
  lapply(group_levels, function(g) {
    Matrix::rowMeans(
      peak_mat[, groups == g, drop = FALSE]
    )
  })
)

colnames(avg_mat) <- group_levels


# -----------------------------------------------------------------------------
# 6. Build column metadata from cluster_stage labels
# -----------------------------------------------------------------------------

col_df <- data.frame(
  sample = colnames(avg_mat),
  stringsAsFactors = FALSE
)

# Parse stage and cluster from the combined label
col_df$stage   <- sub(".*_(3 dpf|adult)$", "\\1", col_df$sample)
col_df$cluster <- sub("_(3 dpf|adult)$",   "",    col_df$sample)


# -----------------------------------------------------------------------------
# 7. Order columns: 3 dpf first, then adult; clusters in default factor order
# Note: brain clusters are ordered alphabetically by default. To enforce a
#       custom order, define a cluster_order vector and set:
#       col_df$cluster <- factor(col_df$cluster, levels = cluster_order)
# -----------------------------------------------------------------------------

col_df$cluster <- factor(col_df$cluster)
col_df$stage   <- factor(col_df$stage, levels = c("3 dpf", "adult"))

col_df <- col_df[order(col_df$stage, col_df$cluster), ]


# -----------------------------------------------------------------------------
# 8. Reorder matrix columns to match col_df
# -----------------------------------------------------------------------------

avg_mat <- avg_mat[, col_df$sample, drop = FALSE]


# -----------------------------------------------------------------------------
# 9. Remove constant peaks (zero variance across groups)
# -----------------------------------------------------------------------------

avg_mat <- avg_mat[apply(avg_mat, 1, var) > 0, , drop = FALSE]


# -----------------------------------------------------------------------------
# 10. Row-scale accessibility (z-score per peak)
# -----------------------------------------------------------------------------

heat_mat <- t(scale(t(avg_mat)))
heat_mat[is.na(heat_mat)] <- 0   # Replace NaN from zero-SD rows


# -----------------------------------------------------------------------------
# 11. Order rows by the column of maximum accessibility
# -----------------------------------------------------------------------------

heat_mat <- heat_mat[order(apply(heat_mat, 1, which.max)), ]


# -----------------------------------------------------------------------------
# 12. Define annotation colors
# -----------------------------------------------------------------------------

# Load cluster color palette from file (one color per line, ordered to match
# the sorted cluster factor levels in col_df$cluster)
cluster_names  <- levels(col_df$cluster)
brain_colpal  <- readLines("heatmap.USE.txt")

# Safety check: number of colors must match number of clusters
if (length(brain_colpal) != length(cluster_names)) {
  warning(
    "Color count (", length(brain_colpal), ") does not match ",
    "cluster count (", length(cluster_names), "). ",
    "Check heatmap.USE.txt."
  )
}

cluster_colors <- setNames(brain_colpal, cluster_names)

ann_colors <- list(
  stage   = c("3 dpf" = "#00CDCD", "adult" = "grey80"),
  cluster = cluster_colors
)


# -----------------------------------------------------------------------------
# 13. Build column annotation dataframe
# -----------------------------------------------------------------------------

ann <- data.frame(
  stage   = col_df$stage,
  cluster = col_df$cluster,
  row.names = col_df$sample
)


# -----------------------------------------------------------------------------
# 14. Plot heatmap
# -----------------------------------------------------------------------------

column_width_setting <- 20

p5 <- pheatmap(
  heat_mat,
  cellwidth         = column_width_setting,
  annotation_col    = ann,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  fontsize          = 10,
  color             = colorRampPalette(c("#518BA9", "white", "#BD8D76"))(100),
  breaks            = seq(-2, 2, length.out = 101)
)

p5
