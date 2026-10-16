# Greenslade_2026
This repository contains R scripts for analyses in the manuscript: **Single-cell chromatin profiling reveals dynamic regulatory logic and enhancer elements in brain and retina development**, Jessie E. Greenslade, Hemagowri Veeravenkatasubramanian, Marisa Reed, and Bushra Raj.

All raw sequencing files and processed/annotated objects can be downloaded at GEO accession GSE324435.

## Overview

We provide a temporally resolved single-cell chromatin accessibility resource of ~95,000 zebrafish brain and retina nuclei spanning larval (3 dpf), juvenile (21 dpf), and adult stages. We define 213 discrete chromatin states and uncover widespread, cell type-specific chromatin reorganization across development.

## Repository Structure

```
Greenslade_2026/
├── scripts/
│   ├── Consensus.Peakset.Script.R
│   ├── Figure2/
│   │   ├── Fig2A.B.R
│   │   ├── Fig2C_Global.Heatmap.R
│   │   ├── Fig2_brain.DARs.R
│   │   ├── Fig2_retina.DARs.R
│   │   ├── Fig2F_I.R
│   │   ├── LinkPeaks.R
│   │   ├── FigS3A_retina.heatmap.R
│   │   ├── FigS3B_brain.heatmap.R
│   │   └── FigS3C_GO.retina.R
│   └── Figure3/
│       ├── Fig3_RunChromVar.R
│       ├── Fig3A_C_Motif.enrichment.R
│       └── Fig3D.E_S6_Jaccard.Index.R
├── supplementary/
│   ├── peaks.bed
│   ├── singlecell.csv
│   ├── brain.background_GO.xlsx
│   ├── retina.background_GO.xlsx
│   └── dr10_dc_weight_matrices.txt
└── README.md
```
---
## R Scripts

### Preprocessing
- **Consensus.Peakset.Script.R** — Generates a consensus peak set across 3 dpf, 21 dpf, and adult timepoints. Builds ChromatinAssay objects and performs QC filtering steps for each Signac object.

---
### Figure 2 — Chromatin Accessibility Landscapes

| Script | Description |
|---|---|
| `Fig2A.B.R` | Bar graphs of retina cell type proportions (Fig2A) and example brain/retina cell type proportions (Fig2B) |
| `Fig2C_Global.Heatmap.R` | Global differential accessibility (DA) analysis across timepoints; plots DA peaks in a heatmap (Fig2C) |
| `Fig2_brain.DARs.R` | DA peak analysis for the integrated brain object; outputs CSV of DA peaks across cell types |
| `Fig2_retina.DARs.R` | DA peak analysis for the integrated retina object; outputs CSV of DA peaks across cell types |
| `Fig2F_I.R` | Dotplot analysis (Fig2F–G) and Barplot analysis (Fig2H–I) |
| `LinkPeaks.R` | Peak-to-gene correlation for all DA peaks across clusters |
| `FigS3A_retina.heatmap.R` | Plots accessibility of cluster-specific DA peaks for retina (FigS3A) |
| `FigS3B_brain.heatmap.R` | Plots accessibility of cluster-specific DA peaks for brain (FigS3B) |
| `FigS3C_GO.retina.R` | GO enrichment analysis for retina DA-linked genes; outputs results table and barplot of top GO terms (FigS3C) |

---
### Figure 3 — Transcription Factor Motif Enrichment

| Script | Description |
|---|---|
| `Fig3_RunChromVar.R` | Runs ChromVAR on Signac objects to compute per-cell TF motif deviation scores |
| `Fig3A_C_Motif.enrichment.R` | DA analysis on ChromVAR assay to identify motif enrichment across clusters; plots heatmaps and FeaturePlots; co-embeds 3 dpf scATAC and scRNA objects to identify candidate regulators (Fig3A–C) |
| `Fig3D.E_S6_Jaccard.Index.R` | Compares enriched motif profiles across timepoints using Jaccard Index analysis (Fig3D–E, FigS6) |

---
## Supplementary Files
| File | Description |
|---|---|
| `peaks.bed` | Aggregate peak calls defining the initial chromatin landscape |
| `singlecell.csv` | Per-barcode QC metrics (fragment counts, TSS enrichment, etc.) required for cell filtering |
| `brain.background_GO.xlsx` | Universal gene background used in brain GO enrichment analysis |
| `retina.background_GO.xlsx` | Universal gene background used in retina GO enrichment analysis |
| `dr10_dc_weight_matrices.txt` | Zebrafish (*D. rerio*) TFBS weight matrices converted to PWMs for ChromVAR motif analysis |

---
## Fragment File Indexing
`CoveragePlot()` from Signac requires a tabix-indexed fragment file (`.tbi`).

**Option 1:** Download pre-indexed `.tbi` files from GEO accession 
**GSE324435**.

**Option 2:** Generate `.tbi` index files locally in R using `Rsamtools`:
see `Consensus.Peakset.Script.R` for the indexing workflow.
