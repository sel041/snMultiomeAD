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
h151 <- readRDS("~/h151_RNA_seurat_object.rds")

# clustering
h151 <- NormalizeData(h151)
h151 <- FindVariableFeatures(h151, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(h151)
h151 <- ScaleData(h151, features = all.genes)
h151 <- RunPCA(h151, features = VariableFeatures(object = seu_clusters))
h151 <- FindNeighbors(h151, dims = 1:15)
h151 <- FindClusters(h151, resolution = 0.5)
h151 <- RunUMAP(h151, dims = 1:15)

# subset and subcluster cell types
# after subclustering, check predicted cell types and the expression of cell type-specific markers
# remove clusters that contain mixed cell types and express marker genes from multiple cell types
# microglia
sub <- subset(h151, idents=c("2"))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.15)

# astro
sub <- subset(h151, idents=c(23))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:4)
sub <- FindNeighbors(sub, dims = 1:4)
sub <- FindClusters(sub, resolution = 0.1)

# oligo
sub <- subset(h151, idents=c("0","11","25","29"))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.2)

# endo
sub <- subset(h151, idents=c("24","28"))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.1)

# dg
sub <- subset(h151, idents=c(3,4))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.1)

# ca
sub <- subset(h151, idents=c(1,7,19,27))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:5)
sub <- FindNeighbors(sub, dims = 1:5)
sub <- FindClusters(sub, resolution = 0.1)

# sub
sub <- subset(h151, idents=c(14,22,21,5,6,12,26,17))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:6)
sub <- FindNeighbors(sub, dims = 1:6)
sub <- FindClusters(sub, resolution = 0.1)

# inhibitory neurons
sub <- subset(h151, idents=c(9, 15, 20, 18, 8, 16, 10, 13))
DefaultAssay(sub) <- "RNA"
sub[["RNA"]]$data <- as(sub[["RNA"]]$data, Class = "dgCMatrix")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub)
sub <- RunUMAP(sub, dims = 1:6)
sub <- FindNeighbors(sub, dims = 1:6)
sub <- FindClusters(sub, resolution = 0.25)

# give NA values to the low quality cluster cells

# filter the low quality cluster cells
h151 <- subset(h151, idents="NA", invert=TRUE)
