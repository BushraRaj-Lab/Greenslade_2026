# =============================================================================
# Peak-Gene Link Analysis using Signac LinkPeaks
# 
# Description: Computes peak-gene links across the full object using
#              LinkPeaks with GeneActivity ("RNA" assay) as the expression assay.
#
# NOTE: This script uses a brain object (brain.integrated.rds) but is
#       fully interchangeable with a retina object — simply swap the
#       input .rds file and update output file names accordingly.
# =============================================================================

# Load libraries --------------------------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(stringr)
library(BSgenome.Drerio.UCSC.danRer11)

# =============================================================================
# Step 1 — Load Object
# NOTE: Replace "brain.integrated.rds" with "retina.integrated.rds"
#       (or equivalent) to run on retina data — all steps are identical
# =============================================================================

brain <- readRDS("brain.integrated.rds")

# =============================================================================
# Step 2 — Normalize RNA Assay
# LinkPeaks requires a normalized RNA/GeneActivity assay as expression proxy
# Run on FULL object — do NOT subset peaks before linking
# =============================================================================

DefaultAssay(brain) <- "RNA"

brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
brain <- ScaleData(brain)

# =============================================================================
# Step 3 — Compute Region Statistics
# Required by LinkPeaks to calculate GC content and other peak statistics
# =============================================================================

DefaultAssay(brain) <- "ATAC"
Idents(brain)       <- brain$labels

brain <- RegionStats(
  object = brain,
  assay  = "ATAC",
  genome = BSgenome.Drerio.UCSC.danRer11
)

# =============================================================================
# Step 4 — Compute Peak-Gene Links
# Run on FULL object with ALL peaks
# Filter to peaks of interest AFTER linking (see Step 5)
# =============================================================================

brain <- LinkPeaks(
  object           = brain,
  peak.assay       = "ATAC",
  expression.assay = "RNA",
  distance         = 5e5,
  min.cells        = 10
)

# =============================================================================
# Step 5 — Extract and Filter Links
# =============================================================================

# Extract all links as a data frame -------------------------------------------
all_links <- as.data.frame(Links(brain[["ATAC"]]))

# Filter to peak of interest --------------------------------------------------
# Update peak coordinates as needed
#peak_of_interest <- "chr7..."

#peak_links <- all_links[all_links$peak == peak_of_interest, ]

cat("Total links found:              ", nrow(all_links),  "\n")
#cat("Links found for target peak:    ", nrow(peak_links), "\n")

# =============================================================================
# Step 6 — Save Outputs
# =============================================================================

write.csv(all_links,  "all.links.brain.csv",  row.names = FALSE)
#write.csv(peak_links, "peak.links.brain.csv", row.names = FALSE)

cat("Saved: all.links.brain.csv\n")
#cat("Saved: peak.links.brain.csv\n")