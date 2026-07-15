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
# load two datasets
mad <- readRDS("~/mouseAD_RNA_seurat_object.rds")
h151 <- readRDS("~/h151_RNA_seurat_object_filtered.rds")

# downsample original mouse AD dataset
mad_sub <- subset(x = mad, downsample = 20000)

# process h151 dataset                
h151[["RNA"]] <- split(h151[["RNA"]], f = h151$library)
h151 <- NormalizeData(h151)
h151 <- FindVariableFeatures(h151)
h151 <- ScaleData(h151)
h151 <- RunPCA(h151)
h151 <- JoinLayers(h151)

combined <- merge(mad_sub, y = h151, add.cell.ids = c("mouseAD", "H151"), project = "snRNA")
combined$data_type <- gsub("_.*", "", rownames(combined@meta.data))
combined <- JoinLayers(combined)
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$data_type)
combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- SketchData(object = combined, ncells = 10000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(combined) <- "sketch"
combined <- FindVariableFeatures(combined, verbose = F)
combined <- ScaleData(combined, verbose = F)
combined <- RunPCA(combined, verbose = F)

# integrate the datasets
combined <- IntegrateLayers(combined, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
    dims = 1:15, k.anchor = 20, reference = which(Layers(combined, search = "data") %in% c("data.mouseAD")),
    verbose = F)
# cluster the integrated data
combined <- FindNeighbors(combined, reduction = "integrated.rpca", dims = 1:15)
combined <- FindClusters(combined, resolution = 2)
combined <- RunUMAP(combined, reduction = "integrated.rpca", dims = 1:15, return.model = T, verbose = F)

combined[["sketch"]] <- JoinLayers(combined[["sketch"]])
combined[["sketch"]] <- split(combined[["sketch"]], f = combined$data_type)

combined <- ProjectIntegration(object = combined, 
                                   sketched.assay = "sketch", 
                                   assay = "RNA", 
                                   reduction = "integrated.rpca")
combined <- ProjectData(
  object = combined,
  sketched.assay = "sketch",
  assay = "RNA",
  sketched.reduction = "integrated.rpca.full",
  full.reduction = "integrated.rpca.full",
  dims = 1:15#,
  #refdata = list(subclass.full = "subclass")
)

combined <- RunUMAP(combined, reduction = "integrated.rpca.full", dims = 1:15, reduction.name = "umap.full",
    reduction.key = "UMAP.full")

# transfer cell type annotation
refer <- subset(combined, data_type %in% c("mouseAD"))
query <- subset(combined, data_type %in% c("H151"))
DefaultAssay(refer) <- "RNA"
DefaultAssay(query) <- "RNA"
anchors <- FindTransferAnchors(reference = refer, query = query, dims = 1:15,
    reference.reduction = "integrated.rpca.full")
predictions <- TransferData(anchorset = anchors, refdata = refer$final_cluster, dims = 1:15)
query@meta.data$predicted.id <- predictions$predicted.id
query@meta.data$prediction.score.max <- predictions$prediction.score.max
                
