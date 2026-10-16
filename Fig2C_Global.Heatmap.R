# =============================================================================
# Timepoint-Specific DARs Analysis
# Description: Identifies differentially accessible regions (DARs) across
#              developmental timepoints and generates a heatmap
# =============================================================================

# Load libraries --------------------------------------------------------------
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# =============================================================================
# Step 1 — Subset Integrated Object
# =============================================================================
#load full.integrated.rds
#OR
#load sub.int from Fig2C R script and skip to step 2
Idents(integrated)<- integrated$broad

sub.int <- subset(integrated, idents = c(
  "Brain_Radial Glia",
  "Retina_Progenitors",
  "Brain_Progenitors",
  "Retina_Neurons + RPE",
  "Brain_Neurons",
  "Brain_Glia",
  "Retina_Muller Glia"
))

# Set identity and assay ------------------------------------------------------
Idents(sub.int)    <- sub.int$dataset
DefaultAssay(sub.int) <- "ATAC"

# =============================================================================
# Step 2 — Run FindMarkers for Each Timepoint Comparison
# =============================================================================

# Helper function: run FindMarkers and add peak name as column ----------------
run_da <- function(object, ident1, ident2) {
  df <- FindMarkers(
    object   = object,
    ident.1  = ident1,
    ident.2  = ident2,
    test.use = "wilcox",
    min.pct  = 0.1
  )
  df$peak <- rownames(df)  # Store peak names as column (replaces Column1)
  return(df)
}

# Run all three comparisons ---------------------------------------------------
da_3v21    <- run_da(sub.int, "3 dpf",  "21 dpf")  # 3 dpf vs 21 dpf
da_3vA     <- run_da(sub.int, "3 dpf",  "adult")   # 3 dpf vs adult
da_21vA    <- run_da(sub.int, "21 dpf", "adult")   # 21 dpf vs adult

# Save all results to CSV -----------------------------------------------------
write.csv(da_3v21, "sub.peaks_3v21.csv",   row.names = TRUE)
write.csv(da_3vA,  "sub.peaks_3vadult.csv", row.names = TRUE)
write.csv(da_21vA, "sub.peaks_21vadult.csv", row.names = TRUE)

# =============================================================================
# Step 3 — Filter Upregulated Peaks Per Timepoint
#          Replaces manual Excel sheet editing
#          Thresholds: p_val_adj < 0.05, avg_log2FC > 0 (upregulated in ident.1)
#                      p_val_adj < 0.05, avg_log2FC < 0 (upregulated in ident.2)
# =============================================================================

# Adjust these thresholds as needed
PADJ_THRESH <- 0.05
LFC_THRESH  <- 1.0      # log2FC threshold (peaks with |LFC| > this are kept)

# Filter function -------------------------------------------------------------
filter_upreg <- function(df, direction = "up") {
  if (direction == "up") {
    # Upregulated in ident.1 (positive LFC)
    df |> filter(p_val_adj < PADJ_THRESH & avg_log2FC > LFC_THRESH)
  } else {
    # Upregulated in ident.2 (negative LFC)
    df |> filter(p_val_adj < PADJ_THRESH & avg_log2FC < -LFC_THRESH)
  }
}

# --- 3 dpf upregulated peaks -------------------------------------------------
peaks_3upreg_vs21  <- filter_upreg(da_3v21, "up")   # replaces sheet "3upreg.vs.21d"
peaks_3upreg_vsA   <- filter_upreg(da_3vA,  "up")   # replaces sheet "3upreg.vs.A"

# --- 21 dpf upregulated peaks ------------------------------------------------
peaks_21upreg_vs3  <- filter_upreg(da_3v21, "down") # replaces sheet "21upreg.vs.3d"
peaks_21upreg_vsA  <- filter_upreg(da_21vA, "up")   # replaces sheet "21upreg.vs.A"

# --- Adult upregulated peaks -------------------------------------------------
peaks_Aupreg_vs3   <- filter_upreg(da_3vA,  "down") # replaces sheet "Aupreg.vs.3d"
peaks_Aupreg_vs21  <- filter_upreg(da_21vA, "down") # replaces sheet "Aupreg.vs.21d"

# =============================================================================
# Step 4 — Build Feature Lists (replaces reading Excel sheets)
# =============================================================================

# Unique peaks upregulated at each timepoint ----------------------------------
tp1_features <- unique(c(peaks_3upreg_vs21$peak,  peaks_3upreg_vsA$peak))
tp2_features <- unique(c(peaks_21upreg_vs3$peak,  peaks_21upreg_vsA$peak))
tp3_features <- unique(c(peaks_Aupreg_vs3$peak,   peaks_Aupreg_vs21$peak))

# All DAR features for heatmap ------------------------------------------------
all_features <- unique(c(tp1_features, tp2_features, tp3_features))

cat("Total DAR features for heatmap:", length(all_features), "\n")
cat("  3 dpf-specific:", length(tp1_features), "\n")
cat("  21 dpf-specific:", length(tp2_features), "\n")
cat("  Adult-specific:", length(tp3_features), "\n")

# =============================================================================
# Step 5 — Generate Heatmap
# =============================================================================

# Set identity and factor order -----------------------------------------------
DefaultAssay(sub.int) <- "ATAC"
Idents(sub.int)       <- sub.int$dataset
sub.int$dataset       <- factor(sub.int$dataset, levels = c("3 dpf", "21 dpf", "adult"))

# Compute average accessibility -----------------------------------------------
avg_exp <- AggregateExpression(
  sub.int,
  assays   = "ATAC",
  features = all_features,
  group.by = "dataset",
  slot     = "data"
)$ATAC

# Log transform and row z-score -----------------------------------------------
mat <- as.matrix(avg_exp)
mat <- log1p(mat)
mat <- t(scale(t(mat)))

# Plot heatmap ----------------------------------------------------------------
column_width_setting <- 50

p1 <- pheatmap(
  mat,
  cellwidth    = column_width_setting,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  col = colorRamp2(c(-1.0, 0, 1.0), c("#518BA9", "white", "#BD8D76"))
)

# Save outputs ----------------------------------------------------------------
#setwd("...")

jpeg(
  file   = "../file.name.jpeg",
  width  = 5, height = 6, units = "in", res = 300
)
p1
dev.off()

pdf(
  file   = "../file.name.pdf",
  width  = 5, height = 6
)
p1
dev.off()

# =============================================================================
# Step 6 — Extract Clustered Row Order
# =============================================================================
nrow(p1)
ordered_indices      <- row_order(p1)
clustered_row_names  <- rownames(p1@matrix)[ordered_indices]

cat("Number of rows in heatmap:", length(ordered_indices), "\n")

write.table(
  clustered_row_names,
  file      = "Global.heatmap_row_names.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)