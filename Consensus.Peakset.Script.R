
# Ensure BiocManager is available first
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the core packages, forcing the installation to overwrite
# any locked/corrupted files.
BiocManager::install(c("Seurat", "Signac", "BSgenome.Drerio.UCSC.danRer11"), force = TRUE)

#creating a common peak set
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(dplyr)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# read in peak sets
peaks.3 <- read.table(
  file = "aggr.adult_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.21 <- read.table(
  file = "aggr.21dpf_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.adult <- read.table(
  file = "aggr.3dpf_peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.3 <- makeGRangesFromDataFrame(peaks.3)
gr.21 <- makeGRangesFromDataFrame(peaks.21)
gr.adult <- makeGRangesFromDataFrame(peaks.adult)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.3, gr.21, gr.adult))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#output peaks from combined peaks to tab delimited for HOMER annotation
write.csv(combined.peaks, "combined.peaks.csv")

# load metadata
md.3 <- read.table(
  file = "aggr.adult_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.21 <- read.table(
  file = "aggr.21dpf_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.adult <- read.table(
  file = "aggr.3dpf_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# perform an initial filtering of low count cells
md.3 <- md.3[md.3$passed_filters > 500, ]
md.21 <- md.21[md.21$passed_filters > 500, ]
md.adult <- md.adult[md.adult$passed_filters > 500, ]


# create fragment objects
frags.3 <- CreateFragmentObject(
  path = "aggr.adult_fragments.tsv.gz",
  cells = rownames(md.3)
)

frags.21 <- CreateFragmentObject(
  path = "aggr.21dpf_fragments.tsv.gz",
  cells = rownames(md.21)
)

frags.adult <- CreateFragmentObject(
  path = "aggr.3dpf_fragments.tsv.gz",
  cells = rownames(md.adult)
)


dpf3.counts <- FeatureMatrix(
  fragments = frags.3,
  features = combined.peaks,
  cells = rownames(md.3)
)

dpf21.counts <- FeatureMatrix(
  fragments = frags.21,
  features = combined.peaks,
  cells = rownames(md.21)
)

adult.counts <- FeatureMatrix(
  fragments = frags.adult,
  features = combined.peaks,
  cells = rownames(md.adult)
)

#create objects
dpf3_assay <- CreateChromatinAssay(dpf3.counts, fragments = frags.3)
dpf3 <- CreateSeuratObject(dpf3_assay, assay = "ATAC", meta.data=md.3)
#saveRDS(dpf3, "dpf3.obj.rds")

dpf21_assay <- CreateChromatinAssay(dpf21.counts, fragments = frags.21)
dpf21 <- CreateSeuratObject(dpf21_assay, assay = "ATAC", meta.data=md.21)
#saveRDS(dpf21, "dpf21.obj.rds")

adult_assay <- CreateChromatinAssay(adult.counts, fragments = frags.adult)
adult <- CreateSeuratObject(adult_assay, assay = "ATAC", meta.data=md.adult)
#saveRDS(adult, "adult.obj.rds")

###Create genome annotations
gtf <- rtracklayer::import('Danio_rerio.GRCz11.112.chr.gtf')

Annotation(dpf3) <- gtf
Annotation(dpf21)<- gtf
Annotation(adult)<- gtf

gtf.df <- as.data.frame(gtf)
gtf@seqnames
# add chr to chromosomes
gtf.df$seqnames <- paste("chr", gtf.df$seqnames, sep = "")
gr <- GRanges(seqnames = gtf.df$seqnames,ranges = IRanges(start = gtf.df$start, end = gtf.df$end))
gr@elementMetadata <- gtf@elementMetadata

Annotation(dpf3) <- gr
Annotation(dpf21)<- gr
Annotation(adult)<- gr

# change to UCSC style since the data was mapped to UCSC
#seqlevels(gr) <- paste0('chr', seqlevels(gr))

dpf3@assays$ATAC@annotation
dpf21@assays$ATAC@annotation
adult@assays$ATAC@annotation

# Extract specific columns
selected_columns <- c("transcript_id", 
                      "gene_name", 
                      "gene_id", 
                      "gene_biotype", 
                      "type")

extracted_data <- dpf3@assays$ATAC@annotation[,selected_columns]
extracted_data <- dpf21@assays$ATAC@annotation[,selected_columns]
extracted_data <- adult@assays$ATAC@annotation[,selected_columns]

colnames(extracted_data@elementMetadata) <- c("tx_id", 
                                              "gene_name", 
                                              "gene_id", 
                                              "gene_biotype", 
                                              "type")


Annotation(dpf3) <- extracted_data
Annotation(dpf21) <- extracted_data
Annotation(adult) <- extracted_data

head(Annotation(dpf3))

head(Fragments(dpf3)[[1]])

####QC pre-processing
dpf3 <- NucleosomeSignal(object = dpf3) 
dpf3$nucleosome_group <- ifelse(dpf3$nucleosome_signal > 4, 'NS > 4', 'NS < 4') 
dpf3 <- TSSEnrichment(dpf3, fast = FALSE) 


dpf3$high.tss <- ifelse(dpf3$TSS.enrichment > 2, 'High', 'Low') 
TSSPlot(dpf3, group.by = 'high.tss') + NoLegend() 
dpf3$pct_reads_in_peaks <- dpf3$peak_region_fragments / dpf3$passed_filters * 100 
dpf3$blacklist_ratio <- dpf3$blacklist_region_fragments / dpf3$peak_region_fragments 

VlnPlot(object = dpf3,  
        features = c('peak_region_fragments','pct_reads_in_peaks', 'nucleosome_signal', 'TSS.enrichment'), 
        pt.size = 0.1, ncol = 4) 

dpf3 <- subset(x = dpf3, 
               subset = peak_region_fragments > 2000 & 
                 peak_region_fragments < 40000 & 
                 pct_reads_in_peaks > 60 & 
                 blacklist_ratio < 0.005 & 
                 nucleosome_signal < 4 &  
                 TSS.enrichment > 2 
)
dpf3
#saveRDS(dpf3, "dpf3.qc.rds")

####QC pre-processing 21dpf
dpf21 <- NucleosomeSignal(object = dpf21) 
dpf21$nucleosome_group <- ifelse(dpf21$nucleosome_signal > 4, 'NS > 4', 'NS < 4') 
dpf21 <- TSSEnrichment(dpf21, fast = FALSE) 


dpf21$high.tss <- ifelse(dpf21$TSS.enrichment > 2, 'High', 'Low') 
TSSPlot(dpf21, group.by = 'high.tss') + NoLegend() 
dpf21$pct_reads_in_peaks <- dpf21$peak_region_fragments / dpf21$passed_filters * 100 
dpf21$blacklist_ratio <- dpf21$blacklist_region_fragments / dpf21$peak_region_fragments 

VlnPlot(object = dpf21,  
        features = c('peak_region_fragments','pct_reads_in_peaks', 'nucleosome_signal', 'TSS.enrichment'), 
        pt.size = 0.1, ncol = 4) 

dpf21 <- subset(x = dpf21, 
                subset = peak_region_fragments > 2000 & 
                  peak_region_fragments < 25000 & 
                  pct_reads_in_peaks > 60 & 
                  blacklist_ratio < 0.005 & 
                  nucleosome_signal < 4 &  
                  TSS.enrichment > 2 
)
dpf21
#saveRDS(dpf21, "dpf21.qc.rds")

####QC pre-processing adult
adult <- NucleosomeSignal(object = adult) 
adult$nucleosome_group <- ifelse(adult$nucleosome_signal > 4, 'NS > 4', 'NS < 4') 
adult <- TSSEnrichment(adult, fast = FALSE) 

saveRDS(adult, "adult.TSSEnrichment.rds") 


adult$high.tss <- ifelse(adult$TSS.enrichment > 2, 'High', 'Low') 
TSSPlot(adult, group.by = 'high.tss') + NoLegend() 
adult$pct_reads_in_peaks <- adult$peak_region_fragments / adult$passed_filters * 100 
adult$blacklist_ratio <- adult$blacklist_region_fragments / adult$peak_region_fragments 

VlnPlot(object = adult,  
        features = c('peak_region_fragments','pct_reads_in_peaks', 'nucleosome_signal', 'TSS.enrichment'), 
        pt.size = 0.1, ncol = 4) 

adult <- subset(x = adult, 
                subset = peak_region_fragments > 2000 & 
                  peak_region_fragments < 60000 & 
                  pct_reads_in_peaks > 60 & 
                  blacklist_ratio < 0.005 & 
                  nucleosome_signal < 4 &  
                  TSS.enrichment > 2 
)
adult
#saveRDS(adult, "adult.qc.rds")
