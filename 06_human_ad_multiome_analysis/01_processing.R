set.seed(1234)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(Signac)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(GenomicRanges)
library(Biobase)
library(readxl)
library(ggplot2) 
library(ggbeeswarm) 
library(dplyr)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer) 
library(colorspace) 
library(scales) 
library(viridis)

setwd("~/humanAD")

# create seurat objects
h5_files <- list.files(pattern = "*.h5")
file_names <- basename(h5_files) 
library <- sub("_\\.h5$", "", file_names)

h5_read <- lapply(h5_files, Read10X_h5)
h5_seu <- NULL
for (i in 1:length(h5_files)) {
  h5_seu[[i]] <- CreateSeuratObject(
    counts = h5_read[[i]]$`Gene Expression`,
    assay = "RNA"
  )
}

# filter low quality cells based on RNA 
for (i in 1:length(library)) {
  h5_seu[[i]] <- subset(x = h5_seu[[i]],
                           subset = nFeature_RNA > 500 
  ) 
}

# doublet detection using scrublet
for (i in 1:length(library)) {
    doublet <- NULL
    doublet <- scrubDoublets(as.matrix(h5_seu[[i]][["RNA"]]$counts), # previous seurat version: as.matrix(h5_seu[[i]]@assays$RNA@counts),
                             expected_doublet_rate = 0.06,
                             verbose = TRUE)
  h5_seu[[i]]$scrubDoublets <- doublet$scrubDoublets
  h5_seu[[i]]$doublet_scores_obs <- doublet$doublet_scores_obs
}                            

# filter doublets
for (i in 1:length(h5_files)) {
  h5_seu[[i]] <- subset(x = h5_seu[[i]],
                           subset = doublet_scores_obs < 0.2 
  )
}

# create signac objects
h5_frags <- list.files(path = "~/outs", pattern = "_fragments.tsv.gz$", full.names = TRUE)

signac_obj <- NULL
for (i in 1:length(library)) {
  counts <- Read10X_h5(h5_files[i])
  chromatin_assay <-  CreateChromatinAssay(
                      counts = counts$Peak,
                      sep = c(":", "-"),
                      genome = "mm10",
                      fragments = h5_frags[i],
                      min.cells = 1)
    signac_obj[[i]] <- CreateSeuratObject(
        counts = chromatin_assay,
        assay = 'peaks',
        project = 'ATAC')
    signac_obj[[i]]@meta.data$cell_barcode <- rownames(signac_obj[[i]]@meta.data)
    signac_obj[[i]]@meta.data$library_id <- library[i]
}

# identifying a common peak set
combined.peaks <- UnifyPeaks(object.list = signac_obj, mode = "reduce")

# create fragment objects
frags <- NULL
for (i in 1:length(library)) {
    frags[[i]] <- CreateFragmentObject(
        path = h5_frags[i],
        cells = colnames(signac_obj[[i]])
    )
}

# quantify peaks in each dataset
signac_obj.counts <- NULL
for (i in 1:length(library)) {
    signac_obj.counts[[i]] <- FeatureMatrix(
        fragments = frags[[i]],
        features = combined.peaks,
        sep = c(":", "-"),
        cells = colnames(signac_obj[[i]])
    )
}

# merge objects
for (i in 1:length(library)) {
    signac_obj[[i]][['peaks']] <- CreateAssayObject(counts = signac_obj.counts[[i]])
}

had_signac_obj <- merge(signac_obj[[1]], y = signac_obj[2:length(signac_obj)], 
                   add.cell.ids = c(library), project = "humanAD")

# filter low quality cells based on # of ATAC peaks
had_signac_obj <- subset(x = had_signac_obj, subset = nCount_peaks > 1000)

# subset cells passed qc from both RNA and ATAC data
cells_passed <- intersect(rownames(had_seurat_obj@meta.data), rownames(had_signac_obj@meta.data))
had_seurat_obj$qc_pass <- rownames(had_seurat_obj@meta.data) %in% cells_passed
had_signac_obj$qc_pass <- rownames(had_signac_obj@meta.data) %in% cells_passed
had_seurat_obj <- subset(x = had_seurat_obj, subset = qc_pass == TRUE)
had_signac_obj <- subset(x = had_signac_obj, subset = qc_pass == TRUE)
