set.seed(1234)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(GenomicRanges)
library(Biobase)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2) 
library(ggbeeswarm) 
library(dplyr)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer) 
library(colorspace) 
library(scales) 
library(viridis)
library(hdf5r)
library(readxl)

setwd("~/h151")

# create seurat objects
h5_files <- list.files(path = "~/outs", pattern = "\\.h5$", full.names = TRUE)
file_names <- basename(h5_files) 
library <- sub("_\\.h5$", "", file_names)

h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- NULL
for (i in 1:length(h5_files)) {
  h5_seurat[[i]] <- CreateSeuratObject(
    counts = h5_read[[i]]$`Gene Expression`,
    assay = "RNA"
  )
}

# filter low quality cells based on # of RNA features
for (i in 1:length(library)) {
  h5_seurat[[i]] <- subset(x = h5_seurat[[i]],
                           subset = nFeature_RNA > 500 
  ) 
}

# doublet detection using scrublet 
for (i in 1:length(library)) {
    doublet <- NULL
    doublet <- scrubDoublets(as.matrix(h5_seurat[[i]][["RNA"]]$counts), 
                             expected_doublet_rate = 0.06,
                             verbose = TRUE)
    h5_seurat[[i]]$scrubDoublets <- doublet$scrubDoublets
    h5_seurat[[i]]$doublet_scores_obs <- doublet$doublet_scores_obs
    }
# filter doublets
for (i in 1:length(h5_files)) {
  h5_seurat[[i]] <- subset(x = h5_seurat[[i]],
                           subset = doublet_scores_obs < 0.2 
  )
}

h151_seurat_obj <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
             add.cell.ids = c(library), project = "H151")

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
combined.peaks

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

h151_signac_obj <- merge(signac_obj[[1]], y = signac_obj[2:length(signac_obj)], 
                   add.cell.ids = c(library), project = "H151")

# filter low quality cells based on # of ATAC peaks
h151_signac_obj <- subset(x = h151_signac_obj, subset = nCount_peaks > 1000)

# subset cells passed qc from both RNA and ATAC data
cells_passed <- intersect(rownames(h151_seurat_obj@meta.data), rownames(h151_signac_obj@meta.data))
h151_seurat_obj$qc_pass <- rownames(h151_seurat_obj@meta.data) %in% cells_passed
h151_signac_obj$qc_pass <- rownames(h151_signac_obj@meta.data) %in% cells_passed
h151_seurat_obj <- subset(x = h151_seurat_obj, subset = qc_pass == TRUE)
h151_signac_obj <- subset(x = h151_signac_obj, subset = qc_pass == TRUE)
