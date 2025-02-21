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
mad <- readRDS("~/mouseAD_seurat_object_qc_filtered.rds")

# normalizing
mad <- NormalizeData(mad, normalization.method = "LogNormalize", scale.factor = 10000)

# feature selection
mad <- FindVariableFeatures(mad, selection.method = "vst", nfeatures = 4000)
mad <- SketchData(
  object = mad,
  ncells = 50000, 
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# scaling
mad <- FindVariableFeatures(mad)
all.genes <- rownames(mad)
mad <- ScaleData(mad, features = all.genes)

# dimensional reduction and initial clustering
mad <- RunPCA(mad)
mad <- FindNeighbors(mad, dims = 1:15)
mad <- FindClusters(mad, resolution = 1)
mad <- RunUMAP(mad, dims = 1:15, return.model = T)

# extend results to the full datasets
mad <- ProjectData(
  object = mad,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)

# remove clusters with average RNA doublet score > 0.1
df <- aggregate(mad@meta.data$doublet_scores_obs, list(mad@meta.data$cluster_full), FUN=mean) 
mad <- subset(mad, idents = subset(df, x > 0.1)[, 1], invert=T)

# clustering after removing high doublet score clusters
mad <- NormalizeData(mad, normalization.method = "LogNormalize", scale.factor = 10000)

mad <- FindVariableFeatures(mad, selection.method = "vst", nfeatures = 4000)
mad <- SketchData(
  object = mad,
  ncells = 50000, 
  method = "LeverageScore",
  sketched.assay = "sketch"
)

mad <- FindVariableFeatures(mad)
all.genes <- rownames(mad)
mad <- ScaleData(mad, features = all.genes)

mad <- RunPCA(mad)
mad <- FindNeighbors(mad, dims = 1:15)
mad <- FindClusters(mad, resolution = 0.05) 
mad <- RunUMAP(mad, dims = 1:15, return.model = T)

mad <- ProjectData(
  object = mad,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)
