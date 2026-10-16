library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(stringr)
library(pheatmap)
library(tidyr)

retina <- readRDS("retina.integrated.rds")

Idents(retina) <- retina$seurat_clusters

DimPlot(retina, label = T,label.size = 5.5)

##########################################################################
#retina UMAP by dataset Supplemental 
Idents(retina) <- retina$dataset

Idents(retina) <- factor(Idents(retina), levels = c("3 dpf",
                                                    "21 dpf",
                                                    "adult") )

DimPlot(retina, cols = c("#00CDCD", "#FF6A6A","grey80"), alpha = 1)

################################################################################
# Differential Accessibility Peak Summary and DotPlot Generation
#
# Purpose:
#   1. Remove duplicate peak calls within each cluster
#   2. Generate summary statistics by developmental stage
#   3. Calculate cell-type abundances
#   4. Create DotPlots
################################################################################
# ------------------------------------------------------------------------------
# Read Input Data
# ------------------------------------------------------------------------------

merged_df <- read_xlsx("Supplemental.Table3.xlsx", sheet = "retina.DA.peaks")
#or upload "retina.merged_df.csv" output from Fig2.retina.DARs.R script
# ------------------------------------------------------------------------------
# Keep Highest-Enrichment Instance of Each Peak per Cluster
# ------------------------------------------------------------------------------

final_peaks <- merged_df %>%
  group_by(cluster, peak) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Assign Developmental Stage
# ------------------------------------------------------------------------------

final_peaks <- final_peaks %>%
  mutate(
    stage = str_extract(comparison, "3 dpf|21 dpf|adult")
  )

# ------------------------------------------------------------------------------
# Count Peaks per Cluster and Stage
# ------------------------------------------------------------------------------

peaks_long <- final_peaks %>%
  group_by(cluster, stage) %>%
  summarise(
    peak_count = n_distinct(peak),
    .groups = "drop"
  ) %>%
  arrange(cluster, stage)

# Wide-format summary table
peaks_table <- peaks_long %>%
  pivot_wider(
    names_from = stage,
    values_from = peak_count,
    values_fill = 0
  )

print(peaks_table, n = 21)

# ------------------------------------------------------------------------------
# Cell Counts per Cluster and Stage
# ------------------------------------------------------------------------------

cells_table <- table(
  retina$int.labels,
  retina$dataset
)

cells_long <- as.data.frame(cells_table)

colnames(cells_long) <- c(
  "cluster",
  "stage",
  "cell_count"
)

print(cells_table)

# ------------------------------------------------------------------------------
# Create Plotting Data Frame
# ------------------------------------------------------------------------------

name_map <- data.frame(
  cluster = levels(Idents(retina))
)

excluded_clusters <- c(
  "Multiplets?_C13"
)

df_plot <- peaks_long %>%
  left_join(cells_long, by = c("cluster", "stage")) %>%
  left_join(name_map, by = "cluster") %>%
  filter(!cluster %in% excluded_clusters) %>%
  group_by(cluster) %>%
  mutate(
    total_cells_per_ident = sum(cell_count, na.rm = TRUE),
    pct_of_total = (cell_count / total_cells_per_ident) * 100
  ) %>%
  ungroup() %>%
  mutate(
    stage = factor(
      stage,
      levels = c("3 dpf", "21 dpf", "adult")
    )
  )

# Preserve plotting order
cluster_order <- df_plot %>%
  distinct(cluster) %>%
  pull(cluster)

df_plot$cluster <- factor(
  df_plot$cluster,
  levels = cluster_order
)
# ------------------------------------------------------------------------------
# Generate DotPlot
# ------------------------------------------------------------------------------

p2 <- ggplot(
  df_plot,
  aes(x = cluster, y = stage)
) +
  geom_point(
    aes(
      size = peak_count,
      fill = pct_of_total > 10
    ),
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "transparent"
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(1, 15),
    breaks = c(5, 100, 500, 1500, 4000, 7000),
    guide = guide_legend(
      override.aes = list(fill = "black")
    )
  ) +
  theme_bw() +
  labs(
    x = "Cell Type / Cluster",
    y = "Developmental Stage",
    size = "Number of Peaks",
    subtitle =
      "Open circles indicate stages containing <10% of cells for a given identity"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "gray90"
    )
  )

p2

################################################################################
# Cluster Composition by Developmental Stage
#
# Purpose:
# Visualize the proportion of cells contributed by each developmental stage
# within retinal cell clusters.
################################################################################

# ------------------------------------------------------------------------------
# Define developmental stage order
# ------------------------------------------------------------------------------

stage_order <- c(
  "adult",
  "21 dpf",
  "3 dpf"
)

retina$dataset <- factor(
  retina$dataset,
  levels = stage_order
)

# ------------------------------------------------------------------------------
# Define stage colors
# ------------------------------------------------------------------------------

stage_colors <- c(
  "grey80",   # adult
  "#FF6A6A",  # 21 dpf
  "#00CDCD"   # 3 dpf
)

# ------------------------------------------------------------------------------
# Select clusters to display
# ------------------------------------------------------------------------------

# To display all clusters, replace clusters_to_keep with:
# clusters_to_keep <- unique(retina$labels)

clusters_to_keep <- df_plot %>%
  distinct(cluster) %>%
  pull(cluster)

# ------------------------------------------------------------------------------
# Prepare metadata for plotting
# ------------------------------------------------------------------------------

meta_filtered <- retina@meta.data %>%
  filter(int.labels %in% clusters_to_keep) %>%
  mutate(
    labels = factor(
      int.labels,
      levels = cluster_order
    )
  )

# ------------------------------------------------------------------------------
# Generate proportional stacked bar plot
# ------------------------------------------------------------------------------

p_bar <- ggplot(
  meta_filtered,
  aes(
    x = labels,
    fill = dataset
  )
) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    values = stage_colors
  ) +
  labs(
    x = "Cluster",
    y = "Proportion of Cells",
    fill = "Developmental Stage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Display plot
p_bar

################################################################################
# Relationship Between Cell Number and DA Peak Counts
#
# Purpose:
# Evaluate whether the number of differentially accessible (DA) peaks is
# associated with the number of cells present in each cluster-stage combination.
################################################################################

# ------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------

combined <- left_join(
  cells_long,
  peaks_long,
  by = c("cluster", "stage")
)

combined_clean <- combined %>%
  filter(
    !is.na(cell_count),
    !is.na(peak_count),
    cell_count > 0,
    peak_count > 0
  )

# ------------------------------------------------------------------------------
# Spearman Correlation
# Non-parametric test of monotonic association
# ------------------------------------------------------------------------------

spearman_result <- cor.test(
  combined_clean$cell_count,
  combined_clean$peak_count,
  method = "spearman"
)

spearman_result

# ------------------------------------------------------------------------------
# Linear Regression
# Effect of cell number on DA peak count
# ------------------------------------------------------------------------------

model_cellcount <- lm(
  log10(peak_count) ~ log10(cell_count),
  data = combined_clean
)

summary(model_cellcount)

# ------------------------------------------------------------------------------
# Scatter Plot
# ------------------------------------------------------------------------------

p_correlation <- ggplot(
  combined_clean,
  aes(
    x = cell_count,
    y = peak_count
  )
) +
  geom_point(
    aes(color = stage),
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    title = "Relationship Between Cell Number and DA Peak Count",
    x = "Cells per Cluster",
    y = "DA Peaks"
  )

p_correlation

# ------------------------------------------------------------------------------
# Linear Regression Controlling for Developmental Stage
# ------------------------------------------------------------------------------

model_stage_adjusted <- lm(
  log10(peak_count) ~ log10(cell_count) + stage,
  data = combined_clean
)

summary(model_stage_adjusted)

# Test whether adding stage improves model fit
anova(
  model_cellcount,
  model_stage_adjusted
)

# ------------------------------------------------------------------------------
# Effect of Stage Alone
# ------------------------------------------------------------------------------

model_stage_only <- lm(
  log10(peak_count) ~ stage,
  data = combined_clean
)

summary(model_stage_only)

# ------------------------------------------------------------------------------
# Cell Number × Stage Interaction
#
# Tests whether the relationship between cell number and DA peaks differs
# among developmental stages.
# ------------------------------------------------------------------------------

model_interaction <- lm(
  log10(peak_count) ~ log10(cell_count) * stage,
  data = combined_clean
)

summary(model_interaction)

# Test whether interaction significantly improves fit
anova(
  model_stage_adjusted,
  model_interaction
)

################################################################################
#Brain analysis
################################################################################
brain <- readRDS("brain.integrated.rds")

Idents(brain) <- brain$clusters

DimPlot(brain, label = T,label.size = 5.5)

##########################################################################
#brain UMAP by dataset Supplemental 
Idents(brain) <- brain$dataset

Idents(brain) <- factor(Idents(brain), levels = c("3 dpf",
                                                    "21 dpf",
                                                    "adult") )

DimPlot(brain, cols = c("#00CDCD", "#FF6A6A","grey80"), alpha = 1)

################################################################################
# Differential Accessibility Peak Summary and DotPlot Generation
#
# Purpose:
#   1. Remove duplicate peak calls within each cluster
#   2. Generate summary statistics by developmental stage
#   3. Calculate cell-type abundances
#   4. Create DotPlots
################################################################################
# ------------------------------------------------------------------------------
# Read Input Data
# ------------------------------------------------------------------------------

merged_df <- read_xlsx("Supplemental.Table3.xlsx", sheet = "brain.DA.peaks")
#or upload "brain.merged_df.csv" output from Fig2.brain.DARs.R script
# ------------------------------------------------------------------------------
# Keep Highest-Enrichment Instance of Each Peak per Cluster
# ------------------------------------------------------------------------------

final_peaks <- merged_df %>%
  group_by(cluster, peak) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Assign Developmental Stage
# ------------------------------------------------------------------------------

final_peaks <- final_peaks %>%
  mutate(
    stage = str_extract(comparison, "3 dpf|21 dpf|adult")
  )

# ------------------------------------------------------------------------------
# Count Peaks per Cluster and Stage
# ------------------------------------------------------------------------------

peaks_long <- final_peaks %>%
  group_by(cluster, stage) %>%
  summarise(
    peak_count = n_distinct(peak),
    .groups = "drop"
  ) %>%
  arrange(cluster, stage)

# Wide-format summary table
peaks_table <- peaks_long %>%
  pivot_wider(
    names_from = stage,
    values_from = peak_count,
    values_fill = 0
  )

print(peaks_table, n = 21)

# ------------------------------------------------------------------------------
# Cell Counts per Cluster and Stage
# ------------------------------------------------------------------------------

cells_table <- table(
  brain$cluster.anno,
  brain$dataset
)

cells_long <- as.data.frame(cells_table)

colnames(cells_long) <- c(
  "cluster",
  "stage",
  "cell_count"
)

print(cells_table)

# ------------------------------------------------------------------------------
# Create Plotting Data Frame
# ------------------------------------------------------------------------------

name_map <- data.frame(
  cluster = levels(Idents(brain))
)


df_plot <- peaks_long %>%
  left_join(cells_long, by = c("cluster", "stage")) %>%
  left_join(name_map, by = "cluster") %>%
  #filter(!cluster %in% excluded_clusters) %>%
  group_by(cluster) %>%
  mutate(
    total_cells_per_ident = sum(cell_count, na.rm = TRUE),
    pct_of_total = (cell_count / total_cells_per_ident) * 100
  ) %>%
  ungroup() %>%
  mutate(
    stage = factor(
      stage,
      levels = c("3 dpf", "21 dpf", "adult")
    )
  )

# Preserve plotting order
cluster_order <- df_plot %>%
  distinct(cluster) %>%
  pull(cluster)

df_plot$cluster <- factor(
  df_plot$cluster,
  levels = cluster_order
)
# ------------------------------------------------------------------------------
# Generate DotPlot
# ------------------------------------------------------------------------------

p2 <- ggplot(
  df_plot,
  aes(x = cluster, y = stage)
) +
  geom_point(
    aes(
      size = peak_count,
      fill = pct_of_total > 10
    ),
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "transparent"
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(1, 15),
    breaks = c(1, 100, 500, 1500, 4000, 8000, 16000),
    guide = guide_legend(
      override.aes = list(fill = "black")
    )
  ) +
  theme_bw() +
  labs(
    x = "Cell Type / Cluster",
    y = "Developmental Stage",
    size = "Number of Peaks",
    subtitle =
      "Open circles indicate stages containing <10% of cells for a given identity"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "gray90"
    )
  )

p2

################################################################################
# Cluster Composition by Developmental Stage
#
# Purpose:
# Visualize the proportion of cells contributed by each developmental stage
# within brain cell clusters.
################################################################################

# ------------------------------------------------------------------------------
# Define developmental stage order
# ------------------------------------------------------------------------------

stage_order <- c(
  "adult",
  "21 dpf",
  "3 dpf"
)

brain$dataset <- factor(
  brain$dataset,
  levels = stage_order
)

# ------------------------------------------------------------------------------
# Define stage colors
# ------------------------------------------------------------------------------

stage_colors <- c(
  "grey80",   # adult
  "#FF6A6A",  # 21 dpf
  "#00CDCD"   # 3 dpf
)

# ------------------------------------------------------------------------------
# Select clusters to display
# ------------------------------------------------------------------------------

# To display all clusters, replace clusters_to_keep with:
# clusters_to_keep <- unique(brain$labels)

clusters_to_keep <- df_plot %>%
  distinct(cluster) %>%
  pull(cluster)

# ------------------------------------------------------------------------------
# Prepare metadata for plotting
# ------------------------------------------------------------------------------

meta_filtered <- brain@meta.data %>%
  filter(cluster.anno %in% clusters_to_keep) %>%
  mutate(
    labels = factor(
      cluster.anno,
      levels = cluster_order
    )
  )

# ------------------------------------------------------------------------------
# Generate proportional stacked bar plot
# ------------------------------------------------------------------------------

p_bar <- ggplot(
  meta_filtered,
  aes(
    x = labels,
    fill = dataset
  )
) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    values = stage_colors
  ) +
  labs(
    x = "Cluster",
    y = "Proportion of Cells",
    fill = "Developmental Stage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Display plot
p_bar

################################################################################
# Relationship Between Cell Number and DA Peak Counts
#
# Purpose:
# Evaluate whether the number of differentially accessible (DA) peaks is
# associated with the number of cells present in each cluster-stage combination.
################################################################################

# ------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------

combined <- left_join(
  cells_long,
  peaks_long,
  by = c("cluster", "stage")
)

combined_clean <- combined %>%
  filter(
    !is.na(cell_count),
    !is.na(peak_count),
    cell_count > 0,
    peak_count > 0
  )

# ------------------------------------------------------------------------------
# Spearman Correlation
# Non-parametric test of monotonic association
# ------------------------------------------------------------------------------

spearman_result <- cor.test(
  combined_clean$cell_count,
  combined_clean$peak_count,
  method = "spearman"
)

spearman_result

# ------------------------------------------------------------------------------
# Linear Regression
# Effect of cell number on DA peak count
# ------------------------------------------------------------------------------

model_cellcount <- lm(
  log10(peak_count) ~ log10(cell_count),
  data = combined_clean
)

summary(model_cellcount)

# ------------------------------------------------------------------------------
# Scatter Plot
# ------------------------------------------------------------------------------

p_correlation <- ggplot(
  combined_clean,
  aes(
    x = cell_count,
    y = peak_count
  )
) +
  geom_point(
    aes(color = stage),
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    title = "Relationship Between Cell Number and DA Peak Count",
    x = "Cells per Cluster",
    y = "DA Peaks"
  )

p_correlation

# ------------------------------------------------------------------------------
# Linear Regression Controlling for Developmental Stage
# ------------------------------------------------------------------------------

model_stage_adjusted <- lm(
  log10(peak_count) ~ log10(cell_count) + stage,
  data = combined_clean
)

summary(model_stage_adjusted)

# Test whether adding stage improves model fit
anova(
  model_cellcount,
  model_stage_adjusted
)

# ------------------------------------------------------------------------------
# Effect of Stage Alone
# ------------------------------------------------------------------------------

model_stage_only <- lm(
  log10(peak_count) ~ stage,
  data = combined_clean
)

summary(model_stage_only)

# ------------------------------------------------------------------------------
# Cell Number × Stage Interaction
#
# Tests whether the relationship between cell number and DA peaks differs
# among developmental stages.
# ------------------------------------------------------------------------------

model_interaction <- lm(
  log10(peak_count) ~ log10(cell_count) * stage,
  data = combined_clean
)

summary(model_interaction)

# Test whether interaction significantly improves fit
anova(
  model_stage_adjusted,
  model_interaction
)

