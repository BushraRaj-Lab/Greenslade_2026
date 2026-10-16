###############################################################
# 
#   inputs:dr10_dc_weight_matrices.txt 
#   3dpf.rds
#
###############################################################
# 1. PACKAGE INSTALLATION

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "Signac",
  "chromVAR",
  "TFBSTools",
  "motifmatchr",
  "BSgenome.Drerio.UCSC.danRer11",
  "BiocParallel"
)

cran_packages <- c(
  "Seurat"
)

install_if_missing <- function(pkg, bioc = FALSE) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("Installing:", pkg, "\n")
    
    if (bioc) {
      BiocManager::install(pkg, ask = FALSE, update = TRUE)
    } else {
      install.packages(
        pkg,
        repos = "https://cloud.r-project.org"
      )
    }
    
  } else {
    cat(pkg, "already installed.\n")
  }
  
}

cat("Checking Bioconductor packages...\n")

for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

cat("\nChecking CRAN packages...\n")

for (pkg in cran_packages) {
  install_if_missing(pkg)
}


# 2. LOAD LIBRARIES
suppressPackageStartupMessages({
  
  library(Signac)
  library(Seurat)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Drerio.UCSC.danRer11)
  library(BiocParallel)
  
})

# 3: DEFINE THE PATH ---
big_file_path <- "dr10_dc_weight_matrices.txt" 

# 4: READ AND PARSE THE FILE ---
if (file.exists(big_file_path)) {
  raw_lines <- readLines(big_file_path)
  starts <- grep("^NA", raw_lines)
  ends <- grep("^//", raw_lines)
  
  pfm_list_raw <- lapply(seq_along(starts), function(i) {
    motif_lines <- raw_lines[starts[i]:ends[i]]
    motif_name <- trimws(gsub("^NA", "", motif_lines[1]))
    mat_lines <- motif_lines[grepl("^[0-9]", motif_lines)]
    if(length(mat_lines) == 0) return(NULL)
    mat_data <- read.table(text = mat_lines, header = FALSE)
    final_mat <- t(as.matrix(mat_data[, 2:5]))
    rownames(final_mat) <- c("A", "C", "G", "T")
    
    if (max(final_mat) < 50) { 
      final_mat <- round(final_mat * 100) 
    } else {
      final_mat <- round(final_mat)
    }
    
    PFMatrix(ID = motif_name, name = motif_name, matrixClass = "Unknown", 
             strand = "+", profileMatrix = final_mat)
  })
  
  pfm_list_raw <- Filter(Negate(is.null), pfm_list_raw)
  pfm_list <- do.call(PFMatrixList, pfm_list_raw)
  pwm_list <- toPWM(pfm_list, type = "prob")
  
  message("Success! pwm_list created.")
  
  # --- MEMORY CLEANUP: Only remove if they actually exist ---
  rm(raw_lines, pfm_list_raw, pfm_list)
  gc() 
} else {
  stop("File 'dr10_dc_weight_matrices.txt' not found at the specified path.")
}

#change to .rds obj to run Chromvar on
obj <- readRDS("file.rds")
DefaultAssay(obj) <- "ATAC"

# --- STEP 5: ADD MOTIFS (With Safety Check) ---
if (exists("pwm_list")) {
  obj <- AddMotifs(
    object = obj,
    genome = BSgenome.Drerio.UCSC.danRer11,
    pfm = pwm_list
  )
  # Clean up pwm_list now that it's stored inside the obj object
  rm(pwm_list)
  gc()
} else {
  stop("pwm_list was never created. Check Step 2 for errors.")
}

# --- STEP 6: RUN CHROMVAR ---
library(BiocParallel)
register(SerialParam()) # Forces single-core to prevent RAM spikes

obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Drerio.UCSC.danRer11,
  verbose = TRUE
)

saveRDS(obj, "chromvar.rds")