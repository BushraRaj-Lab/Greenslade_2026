# Greenslade_2026
This repository contains the R scripts used to create and process the Signac objects used for analysis in the manuscript "Single-cell chromatin profiling reveals dynamic regulatory logic and enhancer elements in brain and retina development", Jessie E. Greenslade, Hemagowri Veeravenkatasubramanian, Marisa Reed, and Bushra Raj.

All raw sequencing files and processed/annotated objects can be downloaded at GEO accession GSE324435.

Overview
We provide a temporally resolved single-cell chromatin accessibility atlas of ~95,000 zebrafish brain and retina nuclei spanning larval (3 dpf), juvenile (21 dpf), and adult stages. We define 212 discrete chromatin states and uncover widespread, cell type-specific chromatin reorganization across development.


Consensus Peakset script

This script outlines the generation of a consensus peakset used to construct ChromatinAssay and Signac objects. It details genome annotation, initial Quality Control (QC), and data filtering steps. The workflow utilizes CellRanger aggr outputs; while peak and metadata files are included here, the larger fragment files and H5 files should be retrieved from GEO (Accession: GSE324435).

Supplementary Files

To ensure reproducibility, this repository includes the necessary outputs from the 10x Genomics CellRanger -aggr pipeline:

    peaks.bed: The aggregate peak calls used to define the initial chromatin landscape.

    singlecell.csv: The per-barcode resource containing QC metrics (fragments, TSS enrichment, etc.) required for filtering.
