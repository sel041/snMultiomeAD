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

mad$final_cluster = factor(mad$final_cluster ,
                           levels=c("DG", "CA1", "CA2", "CA3", "SUB_1", "SUB_2", "SUB_3", "SUB-ProS", "Lamp5", "Pvalb", "Sst", "Vip", "Meis2", "Sncg", "Chandelier", "Cajal-Retzius", "NB", "RGL", "Astro", "OPC", "IOL", "NFOL", "MOL", "DAO", "HMG", "DAM", "IFN", "PVM", "Endo","VLMC", "SMC-Peri", "Choroid-plexus"))
# heatmap for cell type marker genes
library(pheatmap)
library(viridis)
Idents(mad) <- "final_cluster"
cluster.averages <- AverageExpression(object = mad, assays = "RNA")
matrix <- as.matrix(as.data.frame(cluster.averages))
colnames(matrix) <- gsub("RNA.","",colnames(matrix))

# cell type marker gene
genes_to_use = c("Slc17a7","Neurod2", # Excitatory neurons
                 "Prox1", "C1ql2", # DG
                 "Lefty1", "Fibcd1", "Wfs1", "Satb2", "Gpr63", # CA1
                 "Crabp1", "Glul", "Scgn",  # CA2
                 "Arhgef26", "Rnf182","Nptx1", # CA2 & CA3
                 "Iyd","Cdh24", # CA3
                 "Tshz2",  "Slc17a6", "Fn1", # SUB  
                 "Gad1","Gad2", # Inhibitory
                 "Lamp5",
                 "Pvalb",
                 "Sst",
                 "Vip",
                 "Meis2",
                 "Sncg",
                 "Vipr2",
                 "Reln", "Nhlh2", "Ndnf", "Trp73", "Lhx5", # CR
                 "Calb2","Dcx", "Eomes", # NB
                 "Sox2", "Vim", "Nes", "Pax6",  # RGL
                 "Gfap", "Aldoc", "Dbx2", # Astro
                 "Slc1a2",
                 "Pdgfra","Cspg4",  # OPC
                 "9630013A20Rik", # IOL
                 "Mag", "Mog", "Mal", "Opalin", # Oligo
                 "Mbp", "Plp1",
                 "Cx3cr1", "Hexb",
                 "Tmem119", "Siglech", "P2ry12", "Csf1r", # HMG
                 "P2ry13", "Trem2",  # DAM
                 "Mrc1", "F13a1", # PVM
                 "Tek", "Pecam1", "Cldn5", "Ocln", "Flt1","Cdh5",  # Endo  
                 "Dcn", "Lum", # VLMC
                 "Slc22a6",
                 "Pdgfrb", # SMC
                 "Kcnj8", # Pericytes
                 "Ttr") # Choroid

matrix.sub <- subset(matrix, rownames(matrix) %in% genes_to_use)
matrix.sub <- matrix.sub[order(match(rownames(matrix.sub), genes_to_use)),]
matrix.sub <- t(matrix.sub)
rownames(matrix.sub) <- gsub("\\.", " ",rownames(matrix.sub))
dim(matrix.sub)

a <- -1
b <- 5
scale <- "column"

# Set larger margin on the left side to ensure all labels are visible
pheatmap(
  mat = matrix.sub, 
  scale = scale, 
  color = viridis(30), 
  breaks = seq(a, b, by = (b - a) / 30), 
  border_color = "black", 
  cluster_cols = FALSE, 
  cluster_rows = FALSE, 
  show_colnames = TRUE, 
  show_rownames = TRUE, 
  angle_col = 90,  
  fontsize = 25,  
  main = "Gene expression of cell type marker genes",
)

