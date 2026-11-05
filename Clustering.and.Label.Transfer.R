# =============================================================================
# snATAC-seq Processing Framework
# Signac / Seurat Pipeline
# =============================================================================
#
# Description:
#   A generic end-to-end framework outlining the steps used to process 
#   the single-nucleus ATAC-seq data using Signac and Seurat. 
#   This script provides a starting point and will require customization 
#   at key steps based on your dataset.
#
# What you will need to customize:
#   - File paths
#   - Clustering resolution 
#   - The broad category mapping — depends on your RNA reference
#   - Cluster annotation maps
#   - Sub-clustering 
#   - Final label names 
#
# Input:
#   - Pre-built Signac object with ATAC assay
#   - Matched scRNA-seq Seurat reference object (for label transfer)
#
# Output:
#   - RDS checkpoints at key steps
#   - CSV tables: DA peaks, gene activity markers, cluster sizes
#   - Annotated Signac object with:
#       * ATAC + RNA (gene activity) assays
#       * UMAP embedding
#       * Cluster labels
#       * Sub-cluster annotations
#       * RNA label transfer predictions
#
# =============================================================================
# 1.  LIBRARIES
# =============================================================================

library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(scCustomize)
library(cowplot)
library(readxl)


# =============================================================================
# 2.  HELPER FUNCTIONS
# =============================================================================

# Performs: TF-IDF → FindTopFeatures → SVD → (DepthCor check) →
# UMAP → FindNeighbors → FindClusters.
# Inspect DepthCor output and exclude dim 1 if it is correlated with
# sequencing depth (this is common and expected in ATAC data).
run_signac_pipeline <- function(obj,
                                dims   = 2:30,
                                res    = 0.8,
                                algo   = 3,
                                cutoff = 'q0') {
  DefaultAssay(obj) <- "ATAC"
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = cutoff)
  obj <- RunSVD(obj)
  print(DepthCor(obj))   # Inspect this — exclude dim 1 if depth-correlated
  obj <- RunUMAP(obj, reduction = 'lsi', dims = dims)
  obj <- FindNeighbors(obj, reduction = 'lsi', dims = dims)
  obj <- FindClusters(obj, algorithm = algo, resolution = res, verbose = FALSE)
  return(obj)
}


# Finds variable features on the reference, computes CCA transfer anchors,
# and transfers Idents() labels from the scRNA reference to the scATAC query object.
run_label_transfer <- function(reference,
                               query,
                               nfeatures = 5000,
                               cca_dims  = 1:40,
                               lsi_dims  = 2:30) {
  DefaultAssay(reference) <- "RNA"
  DefaultAssay(query)     <- "RNA"
  
  reference <- FindVariableFeatures(reference, nfeatures = nfeatures)
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query     = query,
    reduction = 'cca',
    dims      = cca_dims
  )
  
  predicted <- TransferData(
    anchorset        = anchors,
    refdata          = Idents(reference),
    weight.reduction = query[['lsi']],
    dims             = lsi_dims
  )
  
  query <- AddMetaData(query, predicted)
  return(query)
}


# Each row (ATAC cluster) is normalized to sum to 1, showing what fraction
# of cells in that cluster map to each RNA label.
plot_prediction_heatmap <- function(obj,
                                    x_col    = "seurat_clusters",
                                    y_col    = "predicted.id",
                                    low_col  = "#ffffc8",
                                    high_col = "#7d0025") {
  tbl <- table(obj@meta.data[[x_col]], obj@meta.data[[y_col]])
  tbl <- tbl / rowSums(tbl)
  df  <- as.data.frame(tbl)
  
  ggplot(df, aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient(name = "Fraction of cells",
                        low  = low_col, high = high_col) +
    xlab("ATAC cluster") +
    ylab("RNA cell type label") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(angle = 45, size = 7))
}


# Iterates over clusters and replaces Idents with the most frequent
# predicted.id value in that cluster.
assign_majority_id <- function(obj, clust_col = "seurat_clusters") {
  Idents(obj) <- obj@meta.data[[clust_col]]
  for (cl in levels(obj)) {
    cells    <- WhichCells(obj, idents = cl)
    majority <- names(which.max(table(obj$predicted.id[cells])))
    Idents(obj, cells = cells) <- majority
  }
  return(obj)
}


# Removes cells below a prediction score threshold and plots the updated heatmap.
filter_by_prediction_score <- function(obj,
                                       threshold = 0.5,
                                       x_col     = "seurat_clusters") {
  obj_fil <- subset(obj, subset = prediction.score.max > threshold)
  obj_fil@meta.data[[x_col]]  <- factor(obj_fil@meta.data[[x_col]])
  obj_fil$predicted.id        <- factor(obj_fil$predicted.id)
  print(plot_prediction_heatmap(obj_fil, x_col = x_col))
  return(obj_fil)
}


# =============================================================================
# 3.  LOAD DATA
# =============================================================================

# ── Paths — edit these ────────────────────────────────────────────────────────
atac_rds    <- "./data/your_signac_object.rds"   # Pre-built Signac object
rna_ref_rds <- "./data/your_rna_reference.rds"   # Matched scRNA-seq reference
output_dir  <- "./outputs/"


# ── Create output directory ───────────────────────────────────────────────────
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ── Load objects ──────────────────────────────────────────────────────────────
atac <- readRDS(atac_rds)
rna  <- readRDS(rna_ref_rds)

# =============================================================================
# 4.  DIMENSIONALITY REDUCTION & CLUSTERING
# =============================================================================
# Start with a moderate resolution (e.g. 0.8) and inspect the UMAP.
# Re-run with different resolutions until the clustering matches your
# biological expectations.

# ── Parameters — tune to your data ───────────────────────────────────────────
lsi_dims     <- 2:30   # Adjust if DepthCor reveals additional correlated dims
main_res     <- 0.8    # Starting resolution — tune by inspecting UMAP
cluster_algo <- 3      # 3 = Leiden

# ── Run pipeline ─────────────────────────────────────────────────────────────
atac <- run_signac_pipeline(atac,
                            dims = lsi_dims,
                            res  = main_res,
                            algo = cluster_algo)

saveRDS(atac, file.path(output_dir, "atac_clustered.rds"))
DimPlot(atac, label = TRUE)


# =============================================================================
# 5.  GENE ACTIVITY SCORING
# =============================================================================
# # ── Add gene activity as an assay ─────────────────────────────────────────────
# NOTE: This assay is named "RNA" rather than "GeneActivity" for consistency
# with the GEO-deposited objects associated with this pipeline. Gene activity
# scores are NOT RNA expression — they are derived from chromatin accessibility
# at gene bodies and promoters and used as a proxy for gene expression.
gene.activities <- GeneActivity(atac)
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

gene.activities <- GeneActivity(atac)
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

atac <- NormalizeData(atac,
                      assay                = 'RNA',
                      normalization.method = 'LogNormalize',
                      scale.factor         = median(atac$nCount_RNA))

saveRDS(atac, file.path(output_dir, "atac_with_gene_activity.rds"))

# ── Visualize marker genes of interest ───────────────────────────────────────
# Replace with your own marker genes to help identify clusters.
DefaultAssay(atac) <- 'RNA'

FeaturePlot(atac,
            features   = c("gene1", "gene2", "gene3"),  # <-- replace these
            pt.size    = 0.1,
            order      = TRUE,
            min.cutoff = 'q5',
            max.cutoff = 'q95')


# =============================================================================
# 6.  LABEL TRANSFER FROM scRNA-seq REFERENCE
# =============================================================================
# This transfers cell type labels from a matched scRNA-seq object onto
# your ATAC clusters. The quality of transfer depends on how well the
# RNA reference matches your ATAC dataset (same tissue, same timepoint).
#
# Clusters that map cleanly to one RNA label are well-resolved.
# Clusters with mixed signal may need sub-clustering.
# Label transfer parameters
n_hvf    <- 5000   # Variable features for label transfer
cca_dims <- 1:40   # CCA dims for anchor finding — reduce if reference is small

atac <- run_label_transfer(reference = rna,
                           query     = atac,
                           nfeatures = n_hvf,
                           cca_dims  = cca_dims,
                           lsi_dims  = lsi_dims)

saveRDS(atac, file.path(output_dir, "atac_label_transferred.rds"))

# ── Visualize and validate ────────────────────────────────────────────────
DimPlot(atac, group.by = "predicted.id")
plot_prediction_heatmap(atac, x_col = "seurat_clusters")

# Assign each cluster its most common predicted label
atac <- assign_majority_id(atac, clust_col = "seurat_clusters")

# Filter to high-confidence predictions only
pred_score_cutoff <- 0.5   # Raise this if you want stricter filtering
atac.fil <- filter_by_prediction_score(atac, threshold = pred_score_cutoff)


# =============================================================================
# 7.  CLUSTER ANNOTATION
# =============================================================================
# Use the label transfer and gene activity scores of
# known marker genes to assign a biological identity
# to each cluster.
#
# IMPORTANT: Cluster numbers change every time you re-cluster. This map
# must be re-derived from scratch after any re-clustering. Do not reuse
# a map from a previous run.
#
# Workflow:
#   1. Run plot_prediction_heatmap() and note the dominant RNA label
#      for each ATAC cluster
#   2. Visualize marker genes per expected cell type with FeaturePlot()
#   3. Fill in the map below based on steps 1-2
#   4. Rename and save

# ── Fill in after inspecting your data ───────────────────────────────────────
cluster_annotations <- c(
  "0"  = "CellType_A",
  "1"  = "CellType_B",
  "2"  = "CellType_C"
  # ... add one entry per cluster in your data
)

DefaultAssay(atac) <- "ATAC"
Idents(atac)       <- "seurat_clusters"
atac <- RenameIdents(atac, cluster_annotations)
atac$cluster_anno <- Idents(atac)

DimPlot(atac, group.by = "cluster_anno", label = TRUE, repel = TRUE) +
  NoLegend()
saveRDS(atac, file.path(output_dir, "atac_annotated.rds"))

# Initialize sub.cluster column — will be updated in Section 8 if needed
Idents(atac)      <- "cluster_anno"
atac$sub.cluster  <- as.character(Idents(atac))


# =============================================================================
# 8.  SUB-CLUSTERING (OPTIONAL)
# =============================================================================
# Run this section if certain clusters show mixed cell type signals
# in the label transfer heatmap or marker gene plots.
#
# For each cluster you want to sub-cluster:
#   1. Subset the cells belonging to that cluster
#   2. Re-run the Signac pipeline at a lower resolution
#   3. Check DepthCor again on the subset — sometimes additional dims
#      are depth-correlated in a subset even if they were not in the full data
#   4. Inspect markers and label transfer to annotate the sub-clusters
#   5. Write the new labels back into atac$sub.cluster
#
# Repeat for as many clusters as needed. The example below shows the
# pattern for one cluster — duplicate the block for additional clusters.

# ── Example: sub-clustering one mixed cluster ─────────────────────────────────
res_sub    <- 0.2    # Start low — you are expecting only a few sub-types
sub_dims   <- 2:30   # Adjust after checking DepthCor on this subset

Idents(atac) <- "seurat_clusters"
sub1   <- subset(atac, idents = "0")        # <-- replace "0" with your cluster
sub1.b <- run_signac_pipeline(sub1,
                              dims = sub_dims,
                              res  = res_sub)

# Inspect markers to identify sub-cluster identities
DefaultAssay(sub1.b) <- "RNA"
markers.sub1 <- FindAllMarkers(sub1.b, only.pos = TRUE)
write.csv(markers.sub1, file.path(output_dir, "markers_sub1.csv"))

# Fill in after inspecting markers
sub1_annotations <- c(
  "0" = "SubType_A",
  "1" = "SubType_B"
  # ... add one entry per sub-cluster
)
sub1.b <- RenameIdents(sub1.b, sub1_annotations)
saveRDS(sub1.b, file.path(output_dir, "sub1.rds"))

# Write sub-cluster labels back into the full object
atac$sub.cluster[Cells(sub1.b)] <- as.character(Idents(sub1.b))

# ── Visualize updated sub-cluster labels ─────────────────────────────────────
DimPlot(atac, group.by = "sub.cluster", label = TRUE, repel = TRUE) +
  NoLegend()


# =============================================================================
# 9. FINAL CLUSTER LABELS
# =============================================================================
# Once all annotation and sub-clustering is complete, create a final
# metadata column with clean, consistent cluster names.

# ── Map sub.cluster labels to final names ────────────────────────────────────

final_labels <- c(
  "CellType_A"  = "C0_CellType_A",
  "CellType_B"  = "C1_CellType_B",
  "SubType_A"   = "C2_SubType_A",
  "SubType_B"   = "C3_SubType_B"
  # ... add one entry per final cluster
)

Idents(atac) <- "sub.cluster"
atac <- RenameIdents(atac, final_labels)
atac$Cluster.Labels <- Idents(atac)

# ── Numeric-only labels for programmatic use ─────────────────────────────────
numeric_labels <- setNames(
  sub("_.*", "", as.character(final_labels)),
  final_labels
)
Idents(atac) <- "Cluster.Labels"
atac         <- RenameIdents(atac, numeric_labels)
atac$Clusters <- Idents(atac)

# ── Sort cluster levels numerically ──────────────────────────────────────────
ident_levels  <- levels(atac)
numeric_vals  <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ident_levels)))
sorted_levels <- ident_levels[order(numeric_vals)]
Idents(atac)  <- factor(Idents(atac), levels = sorted_levels)
atac$Clusters <- Idents(atac)

# ── Save final object and cluster size table ──────────────────────────────────
cluster_sizes <- as.data.frame(table(Idents(atac)))
colnames(cluster_sizes) <- c("Cluster", "N_Cells")
write.csv(cluster_sizes,
          file.path(output_dir, "cluster_sizes.csv"),
          row.names = FALSE)

saveRDS(atac, file.path(output_dir, "atac_final.rds"))

DimPlot(atac, label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()


# =============================================================================
# 10. DIFFERENTIAL ACCESSIBILITY & GENE ACTIVITY MARKERS
# =============================================================================

# ──  Differential accessibility peaks (ATAC assay) ───────────────────────
DefaultAssay(atac) <- "ATAC"
Idents(atac)       <- "Cluster.Labels"

da.peaks <- FindAllMarkers(atac, only.pos = TRUE)
da.peaks <- da.peaks %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff), desc(avg_log2FC), .by_group = TRUE)

da.peaks.top1000 <- da.peaks %>%
  group_by(cluster) %>%
  slice_head(n = 1000) %>%
  ungroup()

write.csv(da.peaks,
          file.path(output_dir, "da.peaks_all.csv"),
          row.names = FALSE)
write.csv(da.peaks.top1000,
          file.path(output_dir, "da.peaks_top1000.csv"),
          row.names = FALSE)

# ──  Gene activity markers (RNA assay) ───────────────────────────────────
DefaultAssay(atac) <- "RNA"

gene.markers <- FindAllMarkers(atac, only.pos = TRUE)
gene.markers <- gene.markers %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff), desc(avg_log2FC), .by_group = TRUE)

gene.markers.top1000 <- gene.markers %>%
  group_by(cluster) %>%
  slice_head(n = 1000) %>%
  ungroup()

write.csv(gene.markers,
          file.path(output_dir, "gene_activity_markers_all.csv"),
          row.names = FALSE)
write.csv(gene.markers.top1000,
          file.path(output_dir, "gene_activity_markers_top1000.csv"),
          row.names = FALSE)


# =============================================================================
# END OF SCRIPT
# =============================================================================