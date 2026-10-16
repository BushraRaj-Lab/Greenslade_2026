# =============================================================================
# DAR Accessibility Heatmap
# =============================================================================
# Description: Generates a pseudobulk chromatin accessibility heatmap for
#              differentially accessible regions (DARs) across retinal cell
#              types at 3 dpf and adult stages in zebrafish.
#
# Input:
#   -  retina.integrated.rds 
#   - "retina.merged_df.csv" Output file from Fig2.retina.DARs
#
# Output:
#   - p5: pheatmap object (row-scaled pseudobulk accessibility)
#
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

retina <- readRDS("retina.intgrated.rds")
# -----------------------------------------------------------------------------
# 1. Define exclusions
# -----------------------------------------------------------------------------

# Clusters to exclude entirely (both stages)
exclude_clusters <- c(
  "RGCs/multiplets?_C13"
)

# Specific cluster-stage combinations to exclude
exclude_cluster_stages <- c(
  "Cones_C18_3 dpf",
  "Progenitors_C8_adult",
  "Rods_C16_3 dpf"
)


# -----------------------------------------------------------------------------
# 2. Subset to 3 dpf and adult cells
# -----------------------------------------------------------------------------

retina_plot <- subset(
  retina,
  subset =
    dataset %in% c("3 dpf", "adult") &
    !(int.labels %in% exclude_clusters)
)

# Create a combined cluster_stage label (e.g. "Rods_C5_adult")
retina_plot$cluster_stage <- paste(
  retina_plot$int.labels,
  retina_plot$dataset,
  sep = "_"
)

# Remove specific cluster-stage combinations
retina_plot <- subset(
  retina_plot,
  subset = !(cluster_stage %in% exclude_cluster_stages)
)


# -----------------------------------------------------------------------------
# 3. Select DAR peaks (exclude 21 dpf stage and excluded clusters)
# -----------------------------------------------------------------------------
# Load DAR results
merged_df <- read.csv(
  "retina.merged_df.csv"
)

# Keep highest enrichment peak per cluster
final_peaks <- merged_df %>%
  group_by(cluster, peak) %>%
  slice_max(
    avg_log2FC,
    n = 1,
    with_ties = FALSE
  ) %>% 
  ungroup()%>%
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

DefaultAssay(retina_plot) <- "ATAC"

peak_mat <- GetAssayData(
  retina_plot,
  assay = "ATAC",
  layer = "data"
)[dar_peaks, ]


# -----------------------------------------------------------------------------
# 5. Compute pseudobulk accessibility (mean per cluster_stage group)
# -----------------------------------------------------------------------------

groups       <- retina_plot$cluster_stage
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
col_df$stage <- sub(".*_(3 dpf|adult)$",  "\\1", col_df$sample)
col_df$cluster <- sub("_(3 dpf|adult)$",  "",    col_df$sample)


# -----------------------------------------------------------------------------
# 7. Define cluster display order
# -----------------------------------------------------------------------------

cluster_order <- c(
  "Amacrine Cells_C0",
  "Amacrine Cells_C14",
  "Bipolar Cells_C1",
  "Bipolar Cells_C17",
  "Bipolar Cells_C3",
  "Cones_C15",
  "Cones_C18",
  "Cones_C19",
  "Cones_C2",
  "Cones_C7",
  "Horizontal Cells_C10",
  "Muller Glia_C12",
  "Photoreceptor Precursors_C6",
  "Progenitors_C8",
  "RGCs_C4",
  "RPE_C11",
  "Rods_C16",
  "Rods_C5",
  "Rods_C9"
)


# -----------------------------------------------------------------------------
# 8. Order columns: 3 dpf first, then adult; clusters follow cluster_order
# -----------------------------------------------------------------------------

col_df$cluster <- factor(col_df$cluster, levels = cluster_order)
col_df$stage   <- factor(col_df$stage,   levels = c("3 dpf", "adult"))

col_df <- col_df[order(col_df$stage, col_df$cluster), ]


# -----------------------------------------------------------------------------
# 9. Reorder matrix columns to match col_df
# -----------------------------------------------------------------------------

avg_mat <- avg_mat[, col_df$sample, drop = FALSE]


# -----------------------------------------------------------------------------
# 10. Remove constant peaks (zero variance across groups)
# -----------------------------------------------------------------------------

avg_mat <- avg_mat[apply(avg_mat, 1, var) > 0, , drop = FALSE]


# -----------------------------------------------------------------------------
# 11. Row-scale accessibility (z-score per peak)
# -----------------------------------------------------------------------------

heat_mat <- t(scale(t(avg_mat)))
heat_mat[is.na(heat_mat)] <- 0   # Replace NaN from zero-SD rows


# -----------------------------------------------------------------------------
# 12. Order rows by the column of maximum accessibility
# -----------------------------------------------------------------------------

heat_mat <- heat_mat[order(apply(heat_mat, 1, which.max)), ]


# -----------------------------------------------------------------------------
# 13. Define annotation colors
# -----------------------------------------------------------------------------

celltype_cols <- c(
  "Amacrine Cells_C0"           = "wheat3",
  "Amacrine Cells_C14"          = "tomato3",
  "Bipolar Cells_C1"            = "palegreen3",
  "Bipolar Cells_C17"           = "gold",
  "Bipolar Cells_C3"            = "lightsalmon",
  "Cones_C15"                   = "darkorange",
  "Cones_C18"                   = "pink4",
  "Cones_C19"                   = "palevioletred",
  "Cones_C2"                    = "skyblue1",
  "Cones_C7"                    = "green3",
  "Horizontal Cells_C10"        = "orchid",
  "Muller Glia_C12"             = "mediumpurple2",
  "Photoreceptor Precursors_C6" = "lightseagreen",
  "Progenitors_C8"              = "aquamarine3",
  "RGCs_C4"                     = "goldenrod",
  "RPE_C11"                     = "plum",
  "Rods_C16"                    = "seagreen3",
  "Rods_C5"                     = "rosybrown3",
  "Rods_C9"                     = "khaki3"
)

# Safety check: warn if any cluster is missing a color assignment
missing_colors <- cluster_order[!cluster_order %in% names(celltype_cols)]
if (length(missing_colors) > 0) {
  warning("Missing color entries for: ", paste(missing_colors, collapse = ", "))
}

# Subset colors to only clusters present in the data, preserving order
cluster_colors <- celltype_cols[
  cluster_order[cluster_order %in% levels(col_df$cluster)]
]

ann_colors <- list(
  stage   = c("3 dpf" = "#00CDCD", "adult" = "grey80"),
  cluster = cluster_colors
)


# -----------------------------------------------------------------------------
# 14. Build column annotation dataframe
# -----------------------------------------------------------------------------

ann <- data.frame(
  stage   = col_df$stage,
  cluster = col_df$cluster,
  row.names = col_df$sample
)


# -----------------------------------------------------------------------------
# 15. Plot heatmap
# -----------------------------------------------------------------------------

column_width_setting <- 20

p5 <- pheatmap(
  heat_mat,
  cellwidth        = column_width_setting,
  annotation_col   = ann,
  annotation_colors = ann_colors,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  show_rownames    = FALSE,
  fontsize         = 10,
  color            = colorRampPalette(c("#518BA9", "white", "#BD8D76"))(100),
  breaks           = seq(-2, 2, length.out = 101)
)

p5