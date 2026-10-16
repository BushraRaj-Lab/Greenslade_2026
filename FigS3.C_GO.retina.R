# =============================================================================
# GO ENRICHMENT OF DAR-LINKED GENES & TOP DAR EXPORT
# =============================================================================
# Description: Performs Gene Ontology (GO) enrichment analysis on genes linked
#              to differentially accessible regions (DARs) in zebrafish retina
#              at 3 dpf and adult stages. Produces bar plots per cluster/stage
#              and exports all intermediate and final result tables. 
#              Script interchangeable with brain dataset.  Update background ref
#              with brain.background_GO.xlsx
# Inputs:
#   - retina.background_GO.xlsx          : Background gene universe for GO
#   - merged_df                          : DAR table (in environment) or load
#                                          from "retina.merged_df_09172026.csv"
#   - all.links.retina.csv               : Peak-to-gene link table
#
# Outputs:
#   - retina.DAR_gene_links.csv         : DAR-gene link table
#   - retina.GO_gene_counts.csv         : Gene counts per cluster/stage
#   - retina.GO_gene_lists.csv          : Raw gene lists
#   - retina.genes_entering_GO.csv      : Genes passing GO conversion
#   - retina.GO_results.csv             : Full GO enrichment results
#   - retina.GO_gene_lists_collapsed.csv: Semicolon-collapsed gene lists
#   - retina.GO.3dpf.pdf / .jpeg        : 3 dpf GO enrichment plot
#   - retina.GO.adult.pdf / .jpeg       : Adult GO enrichment plot
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
library(org.Dr.eg.db)


# =============================================================================
# SECTION 3. GO ENRICHMENT OF DAR-LINKED GENES
# =============================================================================
## -----------------------------------------------------------------------------
# 3.1 Build GO Background Universe
# -----------------------------------------------------------------------------

# Load background gene set used as the GO statistical universe.
#
# How this file was generated:
#   A detection rate was calculated for every peak in the ATAC "peaks" assay
#   by dividing the number of cells with a count > 0 by the total number of
#   cells in the retina object. Peaks detected in >= 10% of cells were
#   retained as the high-confidence, reproducibly accessible set.
#   These peaks were then formatted into a three-column BED-style data frame
#   (chrom, start, end) by splitting the peak name on "-" and passed to HOMER
#   for genomic annotation. The resulting HOMER-annotated file, which contains
#   the nearest gene assignments for each peak, was used to extract the unique
#   gene list that serves as the background universe here.
#
# To regenerate:
#   peak_pct <- Matrix::rowSums(
#     GetAssayData(retina, assay = "peaks", layer = "counts") > 0
#   ) / ncol(retina)
#
#   retina_peaks_df <- data.frame(peak = names(peak_pct[peak_pct >= 0.10])) %>%
#     tidyr::separate(peak, into = c("chrom", "start", "end"), sep = "-")
#   # → submit retina_peaks_df to HOMER (annotatePeaks.pl)
#   # → extract unique "Gene Name" column from HOMER output
#   # → save as retina.background_GO.xlsx

ref_genes <- read_excel("retina.background_GO.xlsx")

gene_ref_symbols <- unique(ref_genes$`Gene Name`)

# Convert gene symbols to Entrez IDs using zebrafish OrgDb
converted_universe <- suppressMessages(
  bitr(
    geneID   = gene_ref_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Dr.eg.db
  )
)

entrez_universe_genes <- unique(converted_universe$ENTREZID)

cat("Universe size (Entrez IDs):", length(entrez_universe_genes), "\n")

# -----------------------------------------------------------------------------
# 3.2 Load DARs
# -----------------------------------------------------------------------------

# Uncomment to load from file if merged_df is not already in the environment:
# merged_df <- read.csv("retina.merged_df.csv")

# Deduplicate: keep one row per cluster-peak combination (highest log2FC)
final_peaks <- merged_df %>%
  group_by(cluster, peak) %>%
  slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stage = str_extract(comparison, "3 dpf|21 dpf|adult")
  )


# -----------------------------------------------------------------------------
# 3.3 Add Peak-Gene Links
# -----------------------------------------------------------------------------

# Load peak-to-gene link table (e.g. from Signac LinkPeaks or co-accessibility)
links <- read.csv("all.links.retina.csv")

peaks_link <- final_peaks %>%
  left_join(
    links %>% dplyr::select(peak, gene),
    by = "peak"
  )


# -----------------------------------------------------------------------------
# 3.4 Filter Clusters and Stages
# -----------------------------------------------------------------------------

peaks_link_filtered <- peaks_link %>%
  filter(
    stage != "21 dpf"
  )


# -----------------------------------------------------------------------------
# 3.5 Export DAR-Gene Link Table
# -----------------------------------------------------------------------------

dar_gene_table <- peaks_link_filtered %>%
  filter(!is.na(gene))

write.csv(
  dar_gene_table,
  "retina.DAR_gene_links.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 3.6 Create Per-Cluster/Stage Gene Lists
# -----------------------------------------------------------------------------

gene_lists <- peaks_link_filtered %>%
  filter(!is.na(gene)) %>%
  group_by(cluster, stage) %>%
  summarise(
    genes   = list(unique(gene)),
    n_genes = n_distinct(gene),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 3.7 Export Gene Counts Per Cluster/Stage
# -----------------------------------------------------------------------------

gene_lists %>%
  dplyr::select(cluster, stage, n_genes) %>%
  write.csv("retina.GO_gene_counts.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 3.8 Export Raw Gene Lists (one row per gene)
# -----------------------------------------------------------------------------

gene_lists %>%
  unnest(genes) %>%
  write_csv("retina.GO_gene_lists.csv")


# -----------------------------------------------------------------------------
# 3.9 Gene Symbol → Entrez ID Conversion Function
# -----------------------------------------------------------------------------

convert_genes <- function(gene_symbols,
                          universe_genes = entrez_universe_genes) {
  
  ids <- suppressMessages(
    bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Dr.eg.db
    )
  )
  
  if (is.null(ids) || nrow(ids) == 0) return(NULL)
  
  # Keep only IDs present in the background universe
  ids %>%
    distinct() %>%
    filter(ENTREZID %in% universe_genes)
}


# -----------------------------------------------------------------------------
# 3.10 Apply Gene Conversion Across All Cluster/Stage Groups
# -----------------------------------------------------------------------------

go_input <- gene_lists %>%
  mutate(
    converted = map(genes, ~ convert_genes(.x))
  )


# -----------------------------------------------------------------------------
# 3.11 Export Genes Entering GO Analysis
# -----------------------------------------------------------------------------

# Useful for QC: shows exactly which genes passed symbol → Entrez conversion
genes_entering_go <- go_input %>%
  dplyr::select(cluster, stage, converted) %>%
  unnest(converted)

write.csv(
  genes_entering_go,
  "retina.genes_entering_GO.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 3.12 GO Enrichment Function (Biological Process)
# -----------------------------------------------------------------------------

run_go <- function(converted_ids,
                   universe_genes = entrez_universe_genes) {
  
  # Skip groups with fewer than 5 mapped genes
  if (is.null(converted_ids) || nrow(converted_ids) < 5) return(NULL)
  
  enrichGO(
    gene          = converted_ids$ENTREZID,
    universe      = universe_genes,
    OrgDb         = org.Dr.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",           # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.20,
    readable      = TRUE            # Convert Entrez IDs back to gene symbols
  )
}


# -----------------------------------------------------------------------------
# 3.13 Run GO Enrichment Across All Groups
# -----------------------------------------------------------------------------

go_results <- go_input %>%
  mutate(
    go = map(converted, ~ run_go(.x))
  )


# -----------------------------------------------------------------------------
# 3.14 Extract GO Results Into a Flat Table
# -----------------------------------------------------------------------------

go_df <- go_results %>%
  mutate(
    go_tbl = map(
      go,
      ~ if (is.null(.x)) tibble() else as_tibble(as.data.frame(.x))
    )
  ) %>%
  unnest(go_tbl)


# -----------------------------------------------------------------------------
# 3.15 Export Complete GO Results
# -----------------------------------------------------------------------------

go_df_export <- go_df %>%
  dplyr::select(
    cluster, stage, n_genes,
    ID, Description,
    GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore,
    pvalue, p.adjust, qvalue,
    geneID, Count
  )

write.csv(
  go_df_export,
  "retina.GO_results.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 3.16 Export Gene Lists as Semicolon-Collapsed Strings
# -----------------------------------------------------------------------------

gene_lists_export <- gene_lists %>%
  mutate(
    genes = sapply(genes, paste, collapse = ";")
  )

write.csv(
  gene_lists_export,
  "retina.GO_gene_lists_collapsed.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 3.17 Prepare Data for GO Bar Plots (Top 3 Terms Per Cluster/Stage)
# -----------------------------------------------------------------------------

plot_go <- go_df %>%
  filter(stage %in% c("3 dpf", "adult")) %>%
  group_by(cluster, stage) %>%
  arrange(p.adjust) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    score       = -log10(p.adjust),
    Description = str_wrap(Description, width = 35)
  )


# -----------------------------------------------------------------------------
# Cluster Color Palette
# NOTE: Must match cluster labels exactly — double-check after re-clustering
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
  "Muller Glia_C12"             = "mediumpurple2",   # Note: ensure "_C12" not "_12"
  "Photoreceptor Precursors_C6" = "lightseagreen",
  "Progenitors_C8"              = "aquamarine3",
  "RGCs_C4"                     = "goldenrod",
  "RPE_C11"                     = "plum",
  "Rods_C16"                    = "seagreen3",
  "Rods_C5"                     = "rosybrown3",
  "Rods_C9"                     = "khaki3"
)


# -----------------------------------------------------------------------------
# 3.18 Plot: 3 dpf GO Enrichment
# -----------------------------------------------------------------------------

plot_3dpf <- plot_go %>% filter(stage == "3 dpf")

p_go_3dpf <- ggplot(
  plot_3dpf,
  aes(x = score, y = reorder(Description, score), fill = cluster)
) +
  geom_col() +
  facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = celltype_cols) +
  theme_bw() +
  labs(
    title = "3 dpf GO Enrichment (Biological Process)",
    x     = expression(-log[10](FDR)),
    y     = NULL
  ) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 8),
    legend.position  = "none"
  )

p_go_3dpf

# Save as PDF
pdf("retina.GO.3dpf.pdf", width = 11, height = 15)
print(p_go_3dpf)
dev.off()

# Save as high-resolution JPEG
jpeg("retina.GO.3dpf.jpeg", width = 11, height = 15,
     units = "in", res = 600)
print(p_go_3dpf)
dev.off()


# -----------------------------------------------------------------------------
# 3.19 Plot: Adult GO Enrichment
# -----------------------------------------------------------------------------

plot_adult <- plot_go %>% filter(stage == "adult")

p_go_adult <- ggplot(
  plot_adult,
  aes(x = score, y = reorder(Description, score), fill = cluster)
) +
  geom_col() +
  facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = celltype_cols) +
  theme_bw() +
  labs(
    title = "Adult GO Enrichment (Biological Process)",
    x     = expression(-log[10](FDR)),
    y     = NULL
  ) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 8),
    legend.position  = "none"
  )

p_go_adult

# Save as PDF
pdf("retina.GO.adult.pdf", width = 11, height = 11)
print(p_go_adult)
dev.off()

# Save as high-resolution JPEG
jpeg("retina.GO.adult.jpeg", width = 11, height = 11,
     units = "in", res = 600)
print(p_go_adult)
dev.off()
