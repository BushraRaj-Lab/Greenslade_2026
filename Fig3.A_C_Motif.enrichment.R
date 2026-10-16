# =============================================================================
# Motif Accessibility Heatmap & Co-embedding Analysis
#
# Description:
#   This script performs chromVAR motif accessibility analysis on snATAC-seq
#   data from zebrafish at 3 dpf. It identifies differentially accessible
#   motifs per cell cluster, generates a scaled heatmap, co-embeds scRNA-seq
#   and scATAC-seq data, and produces UMAP feature plots for candidate
#   transcription factors.  This script interchangable with 21 dpf and adult
#   for Heatmap and FeaturePlot sections.scRNA-seq + scATAC-seq Co-embedding
#   requires timepoint and tissue matched scRNA-seq dataset.
#
# Input:
#   - dpf3: A Seurat/Signac object with "chromvar" and "ATAC" assays
#   - rna:  A scRNA-seq dataset
#
# Output:
#   - Motif marker CSV, heatmap PDF, co-embedding RDS, UMAP PDFs/JPEGs
#
# Usage:
#   Set OUTPUT_DIR to your desired output directory before running.
# =============================================================================

library(Seurat)
library(Signac)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)

# -----------------------------------------------------------------------------
# User Configuration
# -----------------------------------------------------------------------------

OUTPUT_DIR <- "Motif.analysis"  # Set your output directory here
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# SECTION 1: Preprocessing — Filter Cell Types
# =============================================================================

Idents(dpf3) <- dpf3$broad
dpf3 <- subset(
  dpf3,
  idents = c("Non Brain Tissue", "Eye_Non Retina"),
  invert = TRUE
)


# =============================================================================
# SECTION 2: Differential Motif Activity (3 dpf)
# =============================================================================

DefaultAssay(dpf3) <- "chromvar"
Idents(dpf3) <- "Cluster.Labels"

# Find differentially active motifs for every cluster
motif_markers_3dpf <- FindAllMarkers(
  object        = dpf3,
  only.pos      = TRUE,
  mean.fxn      = rowMeans,
  fc.name       = "avg_diff",
  logfc.threshold = 0.1,
  test.use      = "wilcox"
)

# Save results
write.csv(
  motif_markers_3dpf,
  file = file.path(OUTPUT_DIR, "motif.markers_3dpf.csv")
)


# =============================================================================
# SECTION 3: Heatmap — Top 25 Motifs per Cluster
# =============================================================================

# Select the top 25 motifs per cluster by average chromVAR deviation
top_motifs_3dpf <- motif_markers_3dpf %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_diff, n = 25) %>%
  pull(gene) %>%
  unique()

# Compute average chromVAR deviation scores per cluster
avg_deviations_3dpf <- AverageExpression(
  object   = dpf3,
  group.by = "Cluster.Labels",
  assays   = "chromvar",
  slot     = "data"
)$chromvar

# Subset to top motifs
heatmap_matrix_3dpf <- avg_deviations_3dpf[top_motifs_3dpf, ]

# --- Scaling: Row-wise Z-score, capped at ±2 ---
matrix_3dpf <- t(scale(t(heatmap_matrix_3dpf)))
matrix_3dpf[matrix_3dpf >  2.0] <-  2.0
matrix_3dpf[matrix_3dpf < -2.0] <- -2.0

# --- Map motif IDs to gene names ---
motif_names_3dpf <- dpf3[["ATAC"]]@motifs@motif.names
rownames(matrix_3dpf) <- motif_names_3dpf[rownames(matrix_3dpf)]

# Retain motif ID for any unmapped entries
unmapped <- is.na(rownames(matrix_3dpf))
rownames(matrix_3dpf)[unmapped] <- names(rownames(matrix_3dpf))[unmapped]

# --- Plot heatmap ---
p <- pheatmap(
  mat           = matrix_3dpf,
  color         = viridis(100, option = "viridis"),
  scale         = "none",
  cluster_rows  = TRUE,
  cluster_cols  = TRUE,
  show_rownames = TRUE,
  border_color  = NA,
  main          = "3 dpf Motif Accessibility (Deviation Score) — Top 25"
)

# Save heatmap to PDF
pdf(
  file   = file.path(OUTPUT_DIR, "Motif.Heatmap.Top25_cluster_cols.FALSE.pdf"),
  width  = 20,
  height = 10
)
print(p)
dev.off()

# --- Save ordered motif list from heatmap clustering ---
row_order        <- p$tree_row$order  
ordered_rownames <- rownames(matrix_3dpf)[row_order]

print(ordered_rownames)
write.csv(
  ordered_rownames,
  file = file.path(OUTPUT_DIR, "rownames.3dpf.motif.top25_cluster_cols.FALSE.csv")
)


# =============================================================================
# SECTION 4: ChromVAR FeaturePlot on dpf3 Object
# =============================================================================

DefaultAssay(dpf3) <- "chromvar"

p_chromvar <- FeaturePlot(
  object    = dpf3,
  features  = "...",
  max.cutoff = "q95",
  min.cutoff = "q5"
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

print(p_chromvar)


# =============================================================================
# SECTION 5: scRNA-seq + scATAC-seq Co-embedding
# =============================================================================
#load 3 dpf scRNA dataset (Raj, 2020) or timepoint/tissue matched scRNA dataset
# --- Prepare ATAC object ---
# Rename the gene activity ("RNA") assay to "ACTIVITY" to avoid collision
dpf3[["ACTIVITY"]] <- dpf3[["RNA"]]
dpf3[["RNA"]]      <- NULL
DefaultAssay(dpf3) <- "ACTIVITY"

# --- Filter both objects to shared cell types ---
Idents(rna)  <- rna$broad
Idents(dpf3) <- dpf3$broad

rna  <- subset(rna,  idents = c("Non Brain Tissue", "Eye Non Retina"), invert = TRUE)
dpf3 <- subset(dpf3, idents = c("Non Brain Tissue", "Eye_Non Retina"), invert = TRUE)

DefaultAssay(rna) <- "RNA"

# --- Find variable features and transfer anchors ---
rna <- FindVariableFeatures(rna, nfeatures = 5000)

transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query     = dpf3,
  reduction = "cca",
  dims      = 1:40
)

# --- Transfer cell type labels to ATAC cells ---
predicted.labels <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = Idents(rna),
  weight.reduction = dpf3[["lsi"]],
  dims             = 2:30
)

dpf3 <- AddMetaData(object = dpf3, metadata = predicted.labels)

# --- Impute RNA expression for ATAC cells ---
rna$dataset  <- "scRNA"
dpf3$dataset <- "scATAC"

genes.use <- VariableFeatures(rna)
refdata   <- GetAssayData(rna, assay = "RNA", layer = "data")[genes.use, ]

imputation   <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = refdata,
  weight.reduction = dpf3[["lsi"]],
  dims             = 2:50
)
dpf3[["RNA"]] <- imputation

# --- Merge and co-embed ---
coembed <- merge(x = rna, y = dpf3)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed,    features = genes.use, verbose  = FALSE)
coembed <- RunUMAP(coembed,   dims     = 1:30)

# --- Visualize co-embedding ---
DimPlot(coembed, group.by = c("dataset", "broad"))
DimPlot(coembed, group.by = "dataset", alpha = 0.3, cols = c("#00CDCD", "#FF6A6A"))

# Save co-embedding object
saveRDS(coembed, file = file.path(OUTPUT_DIR, "coembed.3dpf.rds"))

# --- Plot RNA cell labels ---
rna_cells <- colnames(coembed)[coembed$dataset == "scRNA"]

DimPlot(
  coembed,
  group.by  = "labels",
  cells     = rna_cells,
  label     = FALSE,
  repel     = FALSE,
  reduction = "umap"
) + ggtitle("RNA Cell Labels on Co-embedding")

# --- Plot ATAC cell labels ---
atac_cells <- colnames(coembed)[coembed$dataset == "scATAC"]

DimPlot(
  coembed,
  group.by  = "Cluster.Labels",
  cells     = atac_cells,
  label     = FALSE,
  repel     = FALSE,
  cols      = colpal,  # Ensure `colpal` is defined before running
  reduction = "umap"
) + ggtitle("ATAC Cell Labels on Co-embedding")


# =============================================================================
# SECTION 6: Feature Plots — Candidate Activators
# =============================================================================

save_feature_plot <- function(plot_obj, filename_base, dir = OUTPUT_DIR,
                              width = 11, height = 8) {
  pdf(file.path(dir, paste0(filename_base, ".pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
  jpeg(file.path(dir, paste0(filename_base, ".jpeg")), width = width,
       height = height, units = "in", res = 300)
  print(plot_obj)
  dev.off()
}

# --- RNA expression ---
DefaultAssay(coembed) <- "RNA"

p_rna <- FeaturePlot(
  object     = coembed,
  features   = "gene",
  pt.size    = 0.1,
  order      = TRUE,
  reduction  = "umap",
  min.cutoff = "q5",
  max.cutoff = "q95"
) &
  scale_color_viridis_c(option = "plasma") &
  theme_minimal() &
  theme(panel.grid = element_blank())

print(p_rna)
save_feature_plot(p_rna, "UMAP.rna.gene")

# --- ChromVAR motif accessibility ---
DefaultAssay(coembed) <- "chromvar"
cells_with_motifs <- colnames(coembed[["chromvar"]])

p_motif <- FeaturePlot(
  object     = coembed,
  features   = "gene",
  cells      = cells_with_motifs,
  pt.size    = 0.1,
  order      = TRUE,
  reduction  = "umap",
  min.cutoff = "q5",
  max.cutoff = "q95"
) &
  scale_color_viridis_c(option = "viridis") &
  theme_minimal() &
  theme(panel.grid = element_blank())

print(p_motif)
save_feature_plot(p_motif, "UMAP.motif.gene")

# --- Gene activity ---
DefaultAssay(coembed) <- "ACTIVITY"
cells_with_activity <- colnames(coembed[["ACTIVITY"]])

p_activity <- FeaturePlot(
  object     = coembed,
  features   = "gene",
  cells      = cells_with_activity,
  pt.size    = 0.1,
  order      = TRUE,
  reduction  = "umap",
  min.cutoff = "q5",
  max.cutoff = "q95"
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

print(p_activity)
save_feature_plot(p_activity, "UMAP.gene.activity.gene")