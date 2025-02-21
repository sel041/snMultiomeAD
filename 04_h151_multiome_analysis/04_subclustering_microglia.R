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
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/h151")
# load data
h151 <- readRDS("~/mouseAD_H151_RNA_seurat_object.rds")

# subset microglial cells
Idents(h151) <- "celltype"
mgc <- subset(h151, idents=c("Microgila"))

# subclustering
DefaultAssay(mgc) <- "RNA"
mgc[["RNA"]]$data <- as(mgc[["RNA"]]$data, Class = "dgCMatrix")
mgc <- FindVariableFeatures(mgc)
mgc <- ScaleData(mgc)
mgc <- RunPCA(mgc)
mgc <- RunUMAP(mgc, dims = 1:5)
mgc <- FindNeighbors(mgc, dims = 1:5)
mgc <- FindClusters(mgc, resolution = 0.1)
