set.seed(1234)
library(Signac)
library(Seurat)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

setwd("~/humanAD")
h_sig <- readRDS("~/humanAD_ATAC_signac_object.rds")
hmgc_sig <- subset(h_sig, idents=c("Microgila"))

# get .jaspar file paths
jaspar_dir <- "~/JASPAR2024"
jaspar_files <- list.files(jaspar_dir, pattern = "\\.jaspar$", full.names = TRUE)

# read all files and store PFMatrix objects in a list
pfm_list <- lapply(jaspar_files, function(file) {
  tryCatch({
    readJASPARMatrix(file)
  }, error = function(e) {
    message("Error reading file: ", file, " - ", e$message)
    NULL
  })
})

# remove NULL entries 
pfm_list <- Filter(Negate(is.null), pfm_list)

# combine into a PFMatrixList
pfm_matrix_list <- do.call(c, pfm_list)
pfm_motif_names <- names(pfm_matrix_list)

# set names based on IDs
names(pfm_matrix_list) <- sapply(pfm_list, function(x) ID(x))
pfm <- names(pfm_matrix_list)

# add motif information
hmgc_sig <- AddMotifs(
  object = hmgc_sig,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# computing motif activities
hmgc_sig <- RunChromVAR(
  object = hmgc_sig,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# subclustering based on chromvar motif activities
DefaultAssay(hmgc_sig) <- "chromvar"
hmgc_sig[["chromvar"]]$data <- as(hmgc_sig[["chromvar"]]$data, Class = "dgCMatrix")
hmgc_sig <- FindVariableFeatures(hmgc_sig)
hmgc_sig <- ScaleData(hmgc_sig)
hmgc_sig <- RunPCA(hmgc_sig)
hmgc_sig <- RunUMAP(hmgc_sig, dims = 1:7)
hmgc_sig <- FindNeighbors(hmgc_sig, dims = 1:7)
