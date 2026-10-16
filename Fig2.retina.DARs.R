# Check and install 'future' if missing
if (!requireNamespace("future", quietly = TRUE)) {
  install.packages("future", repos = "https://cloud.r-project.org")
}

# Check and install 'future.apply' if missing
if (!requireNamespace("future.apply", quietly = TRUE)) {
  install.packages("future.apply", repos = "https://cloud.r-project.org")
}

library(future)
library(future.apply)
library(Signac)
library(Seurat)
library(tidyverse)
library(future)

# --- REVISED SPEED CONFIGURATION ---

plan(sequential) 

# 2. Increase the memory limit even further, but be careful
options(future.globals.maxSize = 8000 * 1024^2)

# 3. Explicitly call garbage collection before starting the heavy loop
gc() 

# 1. Setup Variables
obj <- readRDS("retina.integrated.rds") 
timepoint_col <- "dataset"
Idents(obj) <- obj$int.labels
clusters <- levels(Idents(obj))
timepoints <- unique(obj[[timepoint_col, drop = TRUE]])

# Create comparison permutations
comparisons <- expand.grid(ident1 = timepoints, ident2 = timepoints, stringsAsFactors = FALSE) %>%
  filter(ident1 != ident2)

final_list <- list()

comparison_log <- list()


# ... [Upper part of script remains the same] ...

# 2. Iterate through Clusters
for (cluster in clusters) {
  message("--- Processing Cluster: ", cluster, " (", which(clusters == cluster), "/", length(clusters), ") ---")
  
  cluster_sub <- subset(obj, idents = cluster)
  Idents(cluster_sub) <- timepoint_col
  
  # PRE-CHECK: Count cells per dataset in this cluster
  counts <- table(Idents(cluster_sub))
  
  # 3. Iterate through Comparisons (THIS WAS MISSING)
  for (i in 1:nrow(comparisons)) {
    t1 <- comparisons$ident1[i]
    t2 <- comparisons$ident2[i]
    
    # --- ROBUST ERROR FIX ---
    # Get counts safely. Note: counts[t1] from a table can be tricky, 
    # so we check names and convert to numeric.
    c1 <- if(t1 %in% names(counts)) as.numeric(counts[t1]) else 0
    c2 <- if(t2 %in% names(counts)) as.numeric(counts[t2]) else 0
    
    # Log every attempted comparison
    comparison_log[[length(comparison_log) + 1]] <- data.frame(
      cluster = cluster,
      ident1 = t1,
      ident2 = t2,
      cells_ident1 = c1,
      cells_ident2 = c2,
      stringsAsFactors = FALSE
    )
    
    # Ensure values are not NA and meet the 10-cell threshold
    if (!is.na(c1) && !is.na(c2) && c1 >= 10 && c2 >= 10) {
      message("  Comparing: ", t1, " vs ", t2)
      
      markers <- tryCatch({
        FindMarkers(
          object = cluster_sub,
          ident.1 = t1,
          ident.2 = t2,
          test.use = 'wilcox',
          min.pct = 0.1,
          logfc.threshold = 0.5,
          verbose = FALSE
        )
      }, error = function(e) { 
        message("    Error in comparison: ", e$message)
        return(NULL) 
      })
      
      # 4. Filter and Label Results
      if (!is.null(markers) && nrow(markers) > 0) {
        markers_filtered <- markers %>%
          rownames_to_column(var = "peak") %>%
          filter(avg_log2FC > 1.0 & p_val_adj < 0.05) %>%
          mutate(
            cluster = cluster,
            comparison = paste0(t1, "_vs_", t2),
            source_id = paste0("Cluster_", cluster, "_", t1, "_vs_", t2)
          )
        
        final_list[[length(final_list) + 1]] <- markers_filtered
      }
    } else {
      # Log why we skipped
      message("  Skipping ", t1, " vs ", t2, ": Insufficient cells (", c1, " and ", c2, ")")
    }
  } # End of i loop
} # End of cluster loop

comparison_log_df <- bind_rows(comparison_log)

write.csv(
  comparison_log_df,
  "retina.comparison_cell_counts.csv",
  row.names = FALSE
)

# 5. Merge Results
if (length(final_list) > 0) {
  
  merged_df <- bind_rows(final_list)
  
  # Save ALL results before deduplication
  write.csv(
    merged_df,
    "retina.merged_df.csv",
    row.names = FALSE
  )
  
  # Summary from ALL comparisons
  summary_table_raw <- merged_df %>%
    group_by(source_id, cluster, comparison) %>%
    tally(name = "peak_count") %>%
    arrange(desc(peak_count))
  
  write.csv(
    summary_table_raw,
    "retina.summary_raw.csv",
    row.names = FALSE
  )
  
  # Keep strongest instance of each peak
  final_peaks <- merged_df %>%
    group_by(cluster, peak) %>%
    arrange(desc(avg_log2FC)) %>%
    slice(1) %>%
    ungroup()
  
  write.csv(
    final_peaks,
    "retina.final.peaks.csv",
    row.names = FALSE
  )
  
  # Summary after deduplication
  summary_table_final <- final_peaks %>%
    group_by(source_id, cluster, comparison) %>%
    tally(name = "peak_count") %>%
    arrange(desc(peak_count))
  
  write.csv(
    summary_table_final,
    "retina.summary_final.csv",
    row.names = FALSE
  )
  
  print(summary_table_raw)
  print(summary_table_final)
  
} else {
  message("No markers found meeting the criteria.")
}