###############################################################
# Motif Similarity Analysis Across Developmental Timepoints
#
# Assumptions:
# 1. ChromVAR has already been run.
# 2. Three Seurat objects exist in memory:
#      - dpf3
#      - dpf21
#      - adult
# 3. Each object contains:
#      - metadata column: Cluster.Labels
#      - assay: chromvar
# 4. This script creates a metadata column called:
#      - collapsed
###############################################################

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(tibble)
library(proxy)
library(ComplexHeatmap)
library(circlize)

###############################################################
# Helper Functions
###############################################################

collapse_clusters <- function(object, mapping) {
  
  if (!"Cluster.Labels" %in% colnames(object@meta.data)) {
    stop("Metadata column 'Cluster.Labels' not found.")
  }
  
  Idents(object) <- "Cluster.Labels"
  object <- RenameIdents(object, mapping)
  
  object$collapsed <- as.character(Idents(object))
  
  return(object)
}

find_motif_markers <- function(object, timepoint) {
  
  Idents(object) <- object$collapsed
  DefaultAssay(object) <- "chromvar"
  
  markers <- FindAllMarkers(
    object = object,
    only.pos = TRUE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0,
    test.use = "wilcox"
  )
  
  markers$timepoint <- timepoint
  markers$unique_cluster <- paste0(timepoint, "_", markers$cluster)
  
  return(markers)
}

###############################################################
# Cluster Collapse Mappings
###############################################################

dpf3_mapping <- c(
  "C0_Cranial Ganglion" = "Cranial Ganglion",
  "C1_Cranial Ganglion" = "Cranial Ganglion",
  "C2_Cranial Motor Neuron" = "Cranial Ganglion",
  "C3_Diencephalon" = "Diencephalon",
  "C4_Diencephalon" = "Diencephalon",
  "C5_Diencephalon" = "Diencephalon",
  "C6_Diencephalon" = "Diencephalon",
  "C7_Diencephalon" = "Diencephalon",
  "C8_Diencephalon" = "Diencephalon",
  "C9_Dorsal Habenula" = "Habenula",
  "C10_Epidermis" = "Epidermis",
  "C11_Epidermis" = "Epidermis",
  "C12_Epidermis" = "Epidermis",
  "C13_Epithelium" = "Epithelium",
  "C14_Eye_Lens" = "Eye Lens",
  "C15_Granule Cells_Immature" = "Granule Cells",
  "C16_Hindbrain" = "Hindbrain",
  "C17_Hindbrain" = "Hindbrain",
  "C18_Hindbrain" = "Hindbrain",
  "C19_Hindbrain" = "Hindbrain",
  "C20_Hindbrain" = "Hindbrain",
  "C21_Hindbrain" = "Hindbrain",
  "C22_Hindbrain" = "Hindbrain",
  "C23_Hypothalamus" = "Hypothalamus",
  "C24_Hypothalamus" = "Hypothalamus",
  "C25_Iridophores" = "Iridophores",
  "C26_Muscle" = "Muscle",
  "C27_Myeloid" = "Myeloid",
  "C28_Neurons" = "Neurons",
  "C29_Neurons" = "Neurons",
  "C30_Neurons" = "Neurons",
  "C31_Neurons" = "Neurons",
  "C32_Neurons_Glutamatergic" = "Neurons_Glut+",
  "C33_Olfactory Bulb" = "OB",
  "C34_Oligodendrocytes" = "Oligodendrocytes",
  "C35_OPCs" = "OPCs",
  "C36_Optic Tectum" = "Optic Tectum",
  "C37_Optic Tectum" = "Optic Tectum",
  "C38_Otic Vesicle" = "Otic Vesicle",
  "C39_Otic Vesicle" = "Otic Vesicle",
  "C40_Pallium" = "Pallium",
  "C41_Pallium" = "Pallium",
  "C42_Periderm" = "Periderm",
  "C43_Pharyngeal Arch" = "Pharyngeal Arch",
  "C44_Progenitors" = "Progenitors",
  "C45_Progenitors Differentiating" = "Progenitors_Differentiating",
  "C46_Purkinje Neurons" = "Purkinje",
  "C47_Radial Glia" = "Radial Glia",
  "C48_Radial Glia" = "Radial Glia",
  "C49_Retina_Amacrine Cells" = "Amacrine Cells",
  "C50_Retina_Amacrine Cells" = "Amacrine Cells",
  "C51_Retina_Amacrine Cells" = "Amacrine Cells",
  "C52_Retina_Bipolar Cells" = "Bipolar Cells",
  "C53_Retina_Bipolar Cells" = "Bipolar Cells",
  "C54_Retina_Bipolar Cells" = "Bipolar Cells",
  "C55_Retina_Cone Cells" = "Cones",
  "C56_Retina_Cone Cells" = "Cones",
  "C57_Retina_Cone Cells" = "Cones",
  "C58_Retina_Horizontal Cells" = "Horizontal Cells",
  "C59_Retina_Muller Glia" = "Muller Glia",
  "C60_Retina_Photoreceptor Precursors" = "Photoreceptor Precursors",
  "C61_Retina_Photoreceptor Precursors" = "Photoreceptor Precursors",
  "C62_Retina_Progenitors" = "Retina Progenitors",
  "C63_Retina_RGCs" = "RGCs",
  "C64_Retina_RGCs" = "RGCs",
  "C65_Retina_Rod Cells" = "Rods",
  "C66_RPE" = "RPE",
  "C67_RPE" = "RPE",
  "C68_SubPallium" = "SubPallium",
  "C69_Xanthophores" = "Xanthophores"
)

dpf21_mapping <- c(
  "C0 Diencephalon" = "Diencephalon",
  "C1 Diencephalon" = "Diencephalon",
  "C2 Diencephalon" = "Diencephalon",
  "C3 Diencephalon" = "Diencephalon",
  "C4 Diencephalon" = "Diencephalon",
  "C5 Diencephalon" = "Diencephalon",
  "C6 Diencephalon" = "Diencephalon",
  "C7 Diencephalon" = "Diencephalon",
  "C8 Diencephalon" = "Diencephalon",
  "C9 Diencephalon" = "Diencephalon",
  "C10 Dorsal Habenula" = "Habenula",
  "C11 Erythrocytes" = "Erythrocytes",
  "C12 Eye_Lens" = "Eye Lens",
  "C13 Eye_Non-Retina" = "Eye Non-Retina",
  "C14 Eye_Non-Retina" = "Eye Non-Retina",
  "C15 Granule Cells" = "Granule Cells",
  "C16 Granule Cells" = "Granule Cells",
  "C17 Granule Cells" = "Granule Cells",
  "C18 Granule Cells" = "Granule Cells",
  "C19 Hindbrain" = "Hindbrain",
  "C20 Hindbrain_Cranial Ganglion" = "Cranial Ganglion",
  "C21 Hindbrain_Cranial Ganglion" = "Cranial Ganglion",
  "C22 Hindbrain_hox-" = "HB hox-",
  "C23 Hypothalamus" = "Hypothalamus",
  "C24 Hypothalamus" = "Hypothalamus",
  "C25 Hypothalamus" = "Hypothalamus",
  "C26 Immature Neurons" = "IN",
  "C27 Immature Neurons" = "IN",
  "C28 MHB" = "MHB",
  "C29 Microglia" = "Microglia",
  "C30 Neuron pvalb7+_dbh+" = "Neurons_pvalb+",
  "C31 Neuron pvalb7+_dbh+" = "Neurons_pvalb+",
  "C32 Neuron trpc5a+_pvalb6+" = "Neurons_pvalb+",
  "C33 Olfactory Bulb" = "OB",
  "C34 Oligodendrocytes" = "Oligodendrocytes",
  "C35 OPCs" = "OPCs",
  "C36 Optic Tectum" = "Optic Tectum",
  "C37 Optic Tectum" = "Optic Tectum",
  "C38 Optic Tectum" = "Optic Tectum",
  "C39 Optic Tectum" = "Optic Tectum",
  "C40 Pallium" = "Pallium",
  "C41 Pallium" = "Pallium",
  "C42 Pigment Cells" = "Pigment Cells",
  "C43 POA" = "POA",
  "C44 Pretectum/Optic Tectum" = "Pretectum",
  "C45 Progenitors" = "Progenitors",
  "C46 Progenitors" = "Progenitors",
  "C47 Progenitors" = "Progenitors",
  "C48 Progenitors Differentiating (Isthmus)" = "Progenitors Differentiating",
  "C49 Purkinje Neurons" = "Purkinje",
  "C50 Radial Glia" = "Radial Glia",
  "C51 Radial Glia" = "Radial Glia",
  "C52 Retina_Amacrine" = "Amacrine Cells",
  "C53 Retina_Amacrine" = "Amacrine Cells",
  "C54 Retina_Amacrine Cells" = "Amacrine Cells",
  "C55 Retina_Bipolar Cells" = "Bipolar Cells",
  "C56 Retina_Bipolar Cells" = "Bipolar Cells",
  "C57 Retina_Bipolar Cells" = "Bipolar Cells",
  "C58 Retina_Bipolar Cells" = "Bipolar Cells",
  "C59 Retina_Bipolar Cells" = "Bipolar Cells",
  "C60 Retina_Cones" = "Cones",
  "C61 Retina_Cones" = "Cones",
  "C62 Retina_Cones" = "Cones",
  "C63 Retina_Horizontal Cells" = "Horizontal Cells",
  "C64 Retina_Muller Glia" = "Muller Glia",
  "C65 Retina_Photoreceptor Precursors" = "Photoreceptor Precursors",
  "C66 Retina_Photoreceptors" = "Photoreceptors",
  "C67 Retina_Progenitors" = "Retina Progenitors",
  "C68 Retina_RGCs" = "RGCs",
  "C69 Retina_Rods" = "Rods",
  "C70 Retina_RPE" = "RPE",
  "C71 SubPallium" = "SubPallium",
  "C72 TL" = "TL",
  "C73 URL Progenitors" = "URL Progenitors",
  "C74 URL Progenitors" = "URL Progenitors",
  "C75 Ventral Forebrain" = "Ventral FB",
  "C76 Ventral Habenula" = "Habenula"
)

adult_mapping <- c(
  "C0 Cerebellar Neurons" = "Cerebellar Neurons",
  "C1 Diencephalon" = "Diencephalon",
  "C2 Diencephalon" = "Diencephalon",
  "C3 Diencephalon" = "Diencephalon",
  "C4 Diencephalon" = "Diencephalon",
  "C5 Diencephalon" = "Diencephalon",
  "C6 Diencephalon" = "Diencephalon",
  "C7 Diencephalon" = "Diencephalon",
  "C8 Diencephalon" = "Diencephalon",
  "C9 Diencephalon" = "Diencephalon",
  "C10 Epithelium (cornea)" = "Epithelium",
  "C11 Erythrocytes" = "Erythrocytes",
  "C12 Erythrocytes" = "Erythrocytes",
  "C13 Granule Cells" = "Granule Cells",
  "C14 Granule Cells" = "Granule Cells",
  "C15 Granule Cells" = "Granule Cells",
  "C16 Granule Cells" = "Granule Cells",
  "C17 Granule Cells" = "Granule Cells",
  "C18 Granule Cells" = "Granule Cells",
  "C19 Granule Cells" = "Granule Cells",
  "C20 Granule Cells" = "Granule Cells",
  "C21 Granule Cells" = "Granule Cells",
  "C22 Granule Cells" = "Granule Cells",
  "C23 Granule Cells (Immature)" = "Granule Cells",
  "C24 Habenula" = "Habenula",
  "C25 Hindbrain" = "Hindbrain",
  "C26 Hypothalamus" = "Hypothalamus",
  "C27 Hypothalamus" = "Hypothalamus",
  "C28 Hypothalamus" = "Hypothalamus",
  "C29 Hypothalamus" = "Hypothalamus",
  "C30 Immature Neuron" = "IN",
  "C31 Lymphocytes" = "Lymphocytes",
  "C32 Microglia" = "Microglia",
  "C33 Oligodendrocytes" = "Oligodendrocytes",
  "C34 Oligodendrocytes" = "Oligodendrocytes",
  "C35 Oligodendrocytes" = "Oligodendrocytes",
  "C36 Optic Tectum" = "Optic Tectum",
  "C37 Optic Tectum" = "Optic Tectum",
  "C38 Optic Tectum" = "Optic Tectum",
  "C39 Optic Tectum" = "Optic Tectum",
  "C40 Optic Tectum" = "Optic Tectum",
  "C41 Optic Tectum" = "Optic Tectum",
  "C42 Pallium" = "Pallium",
  "C43 Pretectum/Optic Tectum" = "Pretectum",
  "C44 Progenitors" = "Progenitors",
  "C45 Purkinje Cells" = "Purkinje",
  "C46 Radial Glia" = "Radial Glia",
  "C47 Retina_Amacrine Cells" = "Amacrine Cells",
  "C48 Retina_Bipolar Cells" = "Bipolar Cells",
  "C49 Retina_Bipolar Cells" = "Bipolar Cells",
  "C50 Retina_Cones" = "Cones",
  "C51 Retina_Cones" = "Cones",
  "C52 Retina_Cones" = "Cones",
  "C53 Retina_Cones" = "Cones",
  "C54 Retina_Horizontal Cells" = "Horizontal Cells",
  "C55 Retina_Muller Glia" = "Muller Glia",
  "C56 Retina_RGCs" = "RGCs",
  "C57 Retina_Rods" = "Rods",
  "C58 Retina_Rods" = "Rods",
  "C59 Retina_RPE" = "RPE",
  "C60 SubPallium" = "SubPallium",
  "C61 SubPallium" = "SubPallium",
  "C62 SubPallium" = "SubPallium",
  "C63 TL" = "TL",
  "C64 TL" = "TL",
  "C65 URL Progenitors" = "URL Progenitors"
)

###############################################################
# Create collapsed metadata column
###############################################################

dpf3 <- collapse_clusters(dpf3, dpf3_mapping)
saveRDS(dpf3, "dpf3.chromvar.collapsed.rds")

dpf21 <- collapse_clusters(dpf21, dpf21_mapping)
saveRDS(dpf21, "dpf21.chromvar.collapsed.rds")

adult <- collapse_clusters(adult, adult_mapping)
saveRDS(adult, "adult.chromvar.collapsed.rds")

###############################################################
# Motif Marker Discovery
###############################################################

markers_dpf3 <- find_motif_markers(dpf3, "dpf3")
markers_dpf21 <- find_motif_markers(dpf21, "dpf21")
markers_adult <- find_motif_markers(adult, "adult")

write.csv(markers_dpf3, "motif.markers.3dpf.csv", row.names = FALSE)
write.csv(markers_dpf21, "motif.markers.21dpf.csv", row.names = FALSE)
write.csv(markers_adult, "motif.markers.adult.csv", row.names = FALSE)

###############################################################
# Combine marker tables
###############################################################

all_markers <- bind_rows(
  markers_dpf3,
  markers_dpf21,
  markers_adult
)

write.csv(
  all_markers,
  "all.markers.csv",
  row.names = FALSE
)

###############################################################
# Significant motif filtering
###############################################################

sig_motifs <- all_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_diff > 0.8
  ) %>%
  select(
    gene,
    unique_cluster,
    cluster,
    avg_diff
  ) %>%
  distinct()

write.csv(
  sig_motifs,
  "sig.motifs.jaccard.csv",
  row.names = FALSE
)

###############################################################
# Binary motif matrix
###############################################################

binary_matrix <- sig_motifs %>%
  select(-cluster) %>%
  mutate(presence = 1L) %>%
  pivot_wider(
    names_from = gene,
    values_from = presence,
    values_fill = 0L
  ) %>%
  group_by(unique_cluster) %>%
  summarise(across(everything(), max)) %>%
  column_to_rownames("unique_cluster") %>%
  as.matrix()

###############################################################
# Jaccard Similarity
###############################################################

jaccard_sim <- proxy::simil(
  binary_matrix,
  method = "Jaccard"
)

jaccard_mat <- as.matrix(jaccard_sim)
diag(jaccard_mat) <- 1

###############################################################
# Cell types to compare across time
###############################################################

cell_types <- c(
  "Diencephalon",
  "Granule Cells",
  "Habenula",
  "Hindbrain",
  "Hypothalamus",
  "Oligodendrocytes",
  "Optic Tectum",
  "Pallium",
  "Progenitors",
  "Purkinje",
  "Radial Glia",
  "Amacrine Cells",
  "Bipolar Cells",
  "Cones",
  "Horizontal Cells",
  "Muller Glia",
  "RGCs",
  "Rods",
  "RPE",
  "SubPallium"
)

timepoints <- c(
  "dpf3",
  "dpf21",
  "adult"
)

ordered_names <- expand.grid(
  timepoints,
  cell_types
) %>%
  mutate(
    combined = paste(
      Var1,
      Var2,
      sep = "_"
    )
  ) %>%
  pull(combined)

final_order <- ordered_names[
  ordered_names %in% rownames(jaccard_mat)
]

ordered_mat <- jaccard_mat[
  final_order,
  final_order
]

###############################################################
# Global Jaccard Heatmap
###############################################################

heatmap_all <- Heatmap(
  ordered_mat,
  name = "Jaccard",
  col = colorRamp2(
    c(0, 1),
    c("white", "darkblue")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(
    gsub(".*_", "", final_order),
    levels = cell_types
  ),
  column_split = factor(
    gsub(".*_", "", final_order),
    levels = cell_types
  )
)

pdf(
  "Jaccard.All_0.8.pdf",
  width = 20,
  height = 15
)

draw(heatmap_all)

dev.off()

###############################################################
# Example: SubPallium Comparison
###############################################################

target_clusters <- c(
  "dpf3_SubPallium",
  "dpf21_SubPallium",
  "adult_SubPallium"
)

subset_mat <- jaccard_mat[
  target_clusters,
  target_clusters
]

heatmap_subpallium <- Heatmap(
  subset_mat,
  name = "Jaccard Index",
  col = colorRamp2(
    c(0, 1),
    c("white", "darkblue")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

pdf(
  "Jaccard.SubPallium.pdf",
  width = 11,
  height = 8
)

draw(heatmap_subpallium)

dev.off()

###############################################################
# Developmental Motif Conservation and Stage-Specific Motifs
#
# Input:
#   sig_motifs
#
# Required columns:
#   gene
#   cluster
#   unique_cluster
#
# Output:
#   cluster_stage_motif_overlap.csv
###############################################################

library(dplyr)

# Add developmental stage information
motif_df <- sig_motifs %>%
  mutate(
    stage = case_when(
      grepl("^dpf3_", unique_cluster)  ~ "3dpf",
      grepl("^dpf21_", unique_cluster) ~ "21dpf",
      grepl("^adult_", unique_cluster) ~ "adult",
      TRUE ~ NA_character_
    )
  )

# Generate motif lists for each lineage and stage
cluster_lists <- motif_df %>%
  distinct(cluster, stage, gene) %>%
  group_by(cluster, stage) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  ) %>%
  split(.$cluster)

# Compare motif sets across developmental stages
motif_overlap <- lapply(cluster_lists, function(x) {
  
  get_genes <- function(stage_name) {
    
    idx <- which(x$stage == stage_name)
    
    if (length(idx) == 0) {
      character(0)
    } else {
      x$genes[[idx[1]]]
    }
  }
  
  g3  <- get_genes("3dpf")
  g21 <- get_genes("21dpf")
  gad <- get_genes("adult")
  
  list(
    
    # Present at all stages
    shared_all =
      Reduce(intersect, list(g3, g21, gad)),
    
    # Stage-specific motifs
    unique_3dpf =
      setdiff(g3, union(g21, gad)),
    
    unique_21dpf =
      setdiff(g21, union(g3, gad)),
    
    unique_adult =
      setdiff(gad, union(g3, g21)),
    
    # Pairwise shared motifs
    shared_3dpf_21dpf =
      setdiff(intersect(g3, g21), gad),
    
    shared_3dpf_adult =
      setdiff(intersect(g3, gad), g21),
    
    shared_21dpf_adult =
      setdiff(intersect(g21, gad), g3)
  )
})

names(motif_overlap) <- names(cluster_lists)

# Summarize results into a table
motif_overlap_summary <- bind_rows(
  lapply(names(motif_overlap), function(cl) {
    
    data.frame(
      cluster = cl,
      
      n_shared_all =
        length(motif_overlap[[cl]]$shared_all),
      
      n_unique_3dpf =
        length(motif_overlap[[cl]]$unique_3dpf),
      
      n_unique_21dpf =
        length(motif_overlap[[cl]]$unique_21dpf),
      
      n_unique_adult =
        length(motif_overlap[[cl]]$unique_adult),
      
      shared_all =
        paste(motif_overlap[[cl]]$shared_all, collapse = ";"),
      
      unique_3dpf =
        paste(motif_overlap[[cl]]$unique_3dpf, collapse = ";"),
      
      unique_21dpf =
        paste(motif_overlap[[cl]]$unique_21dpf, collapse = ";"),
      
      unique_adult =
        paste(motif_overlap[[cl]]$unique_adult, collapse = ";"),
      
      shared_3dpf_21dpf =
        paste(motif_overlap[[cl]]$shared_3dpf_21dpf, collapse = ";"),
      
      shared_3dpf_adult =
        paste(motif_overlap[[cl]]$shared_3dpf_adult, collapse = ";"),
      
      shared_21dpf_adult =
        paste(motif_overlap[[cl]]$shared_21dpf_adult, collapse = ";")
    )
  })
)

# Save results
write.csv(
  motif_overlap_summary,
  "cluster_stage_motif_overlap.csv",
  row.names = FALSE
)

# View summary in R
motif_overlap_summary