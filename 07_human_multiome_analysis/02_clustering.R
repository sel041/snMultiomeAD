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
had <- readRDS("~/humanAD_RNA_seurat_object.rds")

# initial clustering
had <- NormalizeData(had)
had <- FindVariableFeatures(had)
had <- SketchData(
  object = had,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
DefaultAssay(had) <- "sketch"
had <- FindVariableFeatures(had)
had <- ScaleData(had)
had <- RunPCA(had)
had <- FindNeighbors(had, dims = 1:15)
had <- FindClusters(had, resolution = 1)
had <- RunUMAP(had, dims = 1:15, return.model = T)
had <- ProjectData(
  object = had,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)

# subset and subcluster cell types
# after subclustering, check predicted cell types and the expression of cell type-specific markers
# remove clusters that contain mixed cell types and express marker genes from multiple cell types
# micro
sub <- subset(had, idents=c(3, 7, 10, 28))
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.25)

# oligo
sub <- subset(had, idents=c(2, 11, 29, 0, 1, 26))
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.2)

# excitatory neurons
sub <- subset(had, idents=c(23, 27, 20, 22, 21, 15, 18, 8, 12, 17))
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:6)
sub <- FindNeighbors(sub, dims = 1:6)
sub <- FindClusters(sub, resolution = 0.5)

# inhibitory neurons
sub <- subset(had, idents=c(6, 16, 19, 25, 9, 13)
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:7)
sub <- FindNeighbors(sub, dims = 1:7)
sub <- FindClusters(sub, resolution = 0.2)

# astro
sub <- subset(had, idents=c(4, 14, 5))
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.2)

# endo
sub <- subset(had, idents=c(24))
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.2)
             

