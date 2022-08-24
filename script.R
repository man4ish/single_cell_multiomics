#install.packages("GenomicRanges")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
# read in peak sets
peaks.500 <- read.table(file = "/r_data/pbmc_500/atac_pbmc_500_nextgem_peaks.bed", col.names = c("chr", "start", "end"))
peaks.1k <- read.table(
  file = "/r_data/pbmc_1k/atac_pbmc_1k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.5k <- read.table(
  file = "/r_data/pbmc_5k/atac_pbmc_5k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.10k <- read.table(
  file = "/r_data/pbmc_10k/atac_pbmc_10k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.500 <- makeGRangesFromDataFrame(peaks.500)
gr.1k <- makeGRangesFromDataFrame(peaks.1k)
gr.5k <- makeGRangesFromDataFrame(peaks.5k)
gr.10k <- makeGRangesFromDataFrame(peaks.10k)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.500, gr.1k, gr.5k, gr.10k))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.500 <- read.table(
  file = "/r_data/pbmc_500/atac_pbmc_500_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.1k <- read.table(
  file = "/r_data/pbmc_1k/atac_pbmc_1k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.5k <- read.table(
  file = "/r_data/pbmc_5k/atac_pbmc_5k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.10k <- read.table(
  file = "/r_data/pbmc_10k/atac_pbmc_10k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.500 <- md.500[md.500$passed_filters > 500, ]
md.1k <- md.1k[md.1k$passed_filters > 500, ]
md.5k <- md.5k[md.5k$passed_filters > 500, ]
md.10k <- md.10k[md.10k$passed_filters > 1000, ] # sequenced deeper so set higher cutoff

# create fragment objects
frags.500 <- CreateFragmentObject(
  path = "/r_data/pbmc_500/atac_pbmc_500_nextgem_fragments.tsv.gz",
  cells = rownames(md.500)
)

## Computing hash
frags.1k <- CreateFragmentObject(
  path = "/r_data/pbmc_1k/atac_pbmc_1k_nextgem_fragments.tsv.gz",
  cells = rownames(md.1k)
)
## Computing hash
frags.5k <- CreateFragmentObject(
  path = "/r_data/pbmc_5k/atac_pbmc_5k_nextgem_fragments.tsv.gz",
  cells = rownames(md.5k)
)
## Computing hash
frags.10k <- CreateFragmentObject(
  path = "/r_data/pbmc_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz",
  cells = rownames(md.10k)
)


pbmc500.counts <- FeatureMatrix(
  fragments = frags.500,
  features = combined.peaks,
  cells = rownames(md.500)
)

pbmc1k.counts <- FeatureMatrix(
  fragments = frags.1k,
  features = combined.peaks,
  cells = rownames(md.1k)
)

pbmc5k.counts <- FeatureMatrix(
  fragments = frags.5k,
  features = combined.peaks,
  cells = rownames(md.5k)
)

pbmc10k.counts <- FeatureMatrix(
  fragments = frags.10k,
  features = combined.peaks,
  cells = rownames(md.10k)
)


pbmc500_assay <- CreateChromatinAssay(pbmc500.counts, fragments = frags.500)
pbmc500 <- CreateSeuratObject(pbmc500_assay, assay = "ATAC", meta.data=md.500)

pbmc1k_assay <- CreateChromatinAssay(pbmc1k.counts, fragments = frags.1k)
pbmc1k <- CreateSeuratObject(pbmc1k_assay, assay = "ATAC", meta.data=md.1k)

pbmc5k_assay <- CreateChromatinAssay(pbmc5k.counts, fragments = frags.5k)
pbmc5k <- CreateSeuratObject(pbmc5k_assay, assay = "ATAC", meta.data=md.5k)

pbmc10k_assay <- CreateChromatinAssay(pbmc10k.counts, fragments = frags.10k)
pbmc10k <- CreateSeuratObject(pbmc10k_assay, assay = "ATAC", meta.data=md.10k)

# add information to identify dataset of origin
pbmc500$dataset <- 'pbmc500'
pbmc1k$dataset <- 'pbmc1k'
pbmc5k$dataset <- 'pbmc5k'
pbmc10k$dataset <- 'pbmc10k'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = pbmc500,
  y = list(pbmc1k, pbmc5k, pbmc10k),
  add.cell.ids = c("500", "1k", "5k", "10k")
)
combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
png("/r_data/combined_plot.png")
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
dev.off()
