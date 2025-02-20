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
library(data.table)
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

# merge seurat objects
mad <- merge(seurat_obj[[1]], y = seurat_obj[2:100], 
             add.cell.ids = c(library[1:100]), project = "mouseAD")

# load cell names that passed atac qc (> 1000 reads & > 5 tsse)
pc <- fread("/home/sel041/ps-renlab2/mouseAD/processed_data/atac/qc/atac_qc_passed_cells.csv.gz", header = T)

# keep cells passed quality control from both RNA and ATAC
index <- intersect(rownames(mad@meta.data) , pc$`0`)
mad$qc_pass <- rownames(mad@meta.data) %in% index
mad <- subset(x = mad, subset = qc_pass == TRUE)

# obtain mitochondria reads %
mad[["percent.mt"]] <- PercentageFeatureSet(mad, pattern = "^mt-")
mad <- subset(mad, subset = percent.mt < 10)


