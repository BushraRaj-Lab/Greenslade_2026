# =============================================================================
# Retina and Brain Cell Type Composition Analysis

# Description: Visualizes cell type proportions across developmental stages
#              for zebrafish retina (Figure 2A) and example brain and retina 
#              populations (Figure 2B)
# =============================================================================

# Load libraries --------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# =============================================================================
# Figure 2A — Retina Cell Type Composition (Stacked Bar Chart)
# =============================================================================

# Load data -------------------------------------------------------------------
retina <- read_excel("../aggr/Cluster.sizes_10.12.2026.xlsx", sheet = "retina")

# Reshape: wide → long format -------------------------------------------------
retina_long <- retina |>
  pivot_longer(
    cols      = c(`3 dpf`, `21 dpf`, adult),
    names_to  = "Timepoint",
    values_to = "Count"
  )

# Compute percentages using ALL cells
retina_percent <- retina_long |>
  group_by(Timepoint) |>
  mutate(Percent = Count / sum(Count) * 100)

cell_type_order <- c(
  "Amacrine",
  "Bipolar",
  "Horizontal",
  "RGCs",
  "Cones",
  "Rods",
  "Photoreceptor Precursors",
  "Progenitors",
  "Muller Glia",
  "RPE"
)

# Set factor levels
retina_percent$`Cell Type` <- factor(
  retina_percent$`Cell Type`,
  levels = cell_type_order
)

# Colors
retina_colors <- c(
  "Amacrine"                 = "#FF6A6A",
  "Bipolar"                  = "#FFD700",
  "Horizontal"               = "#99CC00",
  "RGCs"                     = "#FFCCCC",
  "Cones"                    = "#00CDCD",
  "Rods"                     = "#DDA0DD",
  "Photoreceptor Precursors" = "#63B8FF",
  "Progenitors"              = "#FF8C00",
  "Muller Glia"              = "#D3D3D3",
  "RPE"                      = "#DAA520"
)

# Plot
ggplot(
  retina_percent,
  aes(
    x = factor(Timepoint, levels = c("3 dpf", "21 dpf", "adult")),
    y = Percent,
    fill = `Cell Type`
  )
) +
  geom_col(width = 0.4) +
  scale_fill_manual(values = retina_colors) +
  labs(
    x = "Developmental Stage",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
# =============================================================================
# Figure 2B — Brain Cell Type Percentage Change (Line Graph)
# =============================================================================

# Create data frame -----------------------------------------------------------
data <- data.frame(
  CellType = rep(c("Granule Cells", "Progenitors", "Radial Glia",
                   "Rods", "Muller Glia", "Optic Tectum"), each = 3),
  
  Timepoint_Label = rep(c("3 dpf", "21 dpf", "adult"), times = 6),
  
  Percent_Y = c(
    1.56,  9.36,  37.36,   # Granule Cells
    14.34, 9.93,  1.75,    # Progenitors
    5.52,  4.20,  3.21,    # Radial Glia
    1.21,  7.50,  29.40,   # Rods
    1.86,  2.70,  4.98,    # Muller Glia
    7.19,  13.53, 13.67    # Optic Tectum
  )
)

# Set timepoint factor order --------------------------------------------------
data$Timepoint_Label <- factor(data$Timepoint_Label,
                               levels = c("3 dpf", "21 dpf", "adult"))

# Define colors ---------------------------------------------------------------
# Line colors
line_colors <- c(
  "Granule Cells" = "darkgreen",
  "Progenitors"   = "#E69F00",
  "Radial Glia"   = "darkblue",
  "Rods"          = "purple",
  "Muller Glia"   = "brown3",
  "Optic Tectum"  = "deeppink2"
)

# Dot fill colors (slightly darker than line colors)
dot_colors <- c(
  "Granule Cells" = "#006400",
  "Progenitors"   = "#CC8400",
  "Radial Glia"   = "#000080",
  "Rods"          = "#4B0082",
  "Muller Glia"   = "brown4",
  "Optic Tectum"  = "deeppink3"
)

# Plot ------------------------------------------------------------------------
ggplot(data,
       aes(x     = Timepoint_Label,
           y     = Percent_Y,
           group = CellType,
           color = CellType)) +
  geom_line(linewidth = 1) +
  geom_point(
    aes(fill = CellType),
    shape = 21,
    size  = 6
  ) +
  scale_color_manual(values = line_colors) +
  scale_fill_manual(values = dot_colors, guide = "none") +  # hide duplicate legend
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 40)
  ) +
  labs(
    title = "Cell Type Percentage Change Across Developmental Timepoints",
    x     = "Developmental Stage",
    y     = "Percentage (%)",
    color = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )