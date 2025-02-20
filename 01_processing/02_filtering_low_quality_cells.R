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

setwd("~/mouseAD")
# load seurat object
metadata <- read_excel("metadata.xlsx")
library <- metadata$library_id
seurat_obj <- readRDS("~/mouseAD_seurat_object.rds")

# filter low quality cells based on # of RNA features
for (i in 1:length(library)) {
  seurat_obj[[i]] <- subset(x = seurat_obj[[i]],
                           subset = nFeature_RNA > 500 
  ) 
}

# doublet detection using scrublet
# run for each library
doublet <- NULL
doublet <- scrubDoublets(as.matrix(seurat_obj[[i]]@assays$RNA@counts),
                         expected_doublet_rate = 0.06,
                         verbose = TRUE)
saveRDS(doublet, paste0("~/doublet/doublet_",i,".rds"))

# add doublet information
for (i in 1:length(library)) {
  doublet <- readRDS(paste0("~/doublet/doublet_",i,".rds"))
  seurat_obj[[i]]$scrubDoublets <- doublet$scrubDoublets
  seurat_obj[[i]]$doublet_scores_obs <- doublet$doublet_scores_obs
}                            

# filter doublet
for (i in 1:length(library)) {
  seurat_obj[[i]] <- subset(x = h5_seurat[[i]],
                           subset = doublet_scores_obs < 0.2 
  )
}

# load cell name that passed atac qc 



