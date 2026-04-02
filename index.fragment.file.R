
if (!requireNamespace("Rsamtools", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Rsamtools")
}

library(Rsamtools)

# Index the file
Rsamtools::indexTabix("fragments.tsv.gz", format = "bed")

# Create the fragment object (it will now find the .tbi automatically)
frags <- CreateFragmentObject(
  path = "fragments.tsv.gz",
  cells = colnames(your_seurat_object)
)

# Replace the fragment information in your assay
Fragments(your_seurat_object[["ATAC"]]) <- NULL # Clear old path if necessary
Fragments(your_seurat_object[["ATAC"]]) <- frags