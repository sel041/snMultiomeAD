set.seed(1234)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(Signac)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(GenomicRanges)
library(Biobase)
library(readxl)
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

setwd("~/mouseAD/")
mad <- readRDS("~/mouseAD_RNA_seurat_object.rds")

# add gene module score
# genes from Gulen et al., Nature 2023
gene_names <- list(c("Slamf9", "Cd84", "Ly9", "Ctsz", "Cst7", "Cybb", "Cd52", "Cd72", "Gpnmb", "Clec7a", "Cd9", "Itgam", "Itgax", "Lyz2", "Lirb4a", "Lgals3", "Plau", "St14", "Mmp12", "Ccl6", "Cd68", "Cxcl16", "Lgals3bp", "Gpr84", "Tspo", "Apod", "Cd86", "Trem2", "C4b", "Cd74", "Cd274", "Mpeg1"))
mad <- AddModuleScore(
  mad,
  gene_names,
  name = "Activation_Disease_association_score"
)

gene_names <- list(c("B2m", "Tap1", "H2-K1", "H2-D1", "H2-Q7"))
mad <- AddModuleScore(
  mad,
  gene_names,
  name = "MHC_class_I_score"
)

gene_names <- list(c("Stat1", "Sp100", "Ifi204", "Zbp1", "Ifi44", "Isg15", "Oas2", "Oasl2", "Oasl1", "Usp18", "Irf7", "Nlrc5", "Xaf1", "Rnf213", "Ifi27l2a", "Rsad2", "Rtp4", "Ifit2", "Ifit3", "Ifit3b", "Ifit1"))
mad <- AddModuleScore(
  mad,
  gene_names,
  name = "Interferon_related_score"
)

mgc <- readRDS("~/mouseAD_RNA_seurat_object_microglia.rds")
gene_names <- list(c("Slamf9", "Cd84", "Ly9", "Ctsz", "Cst7", "Cybb", "Cd52", "Cd72", "Gpnmb", "Clec7a", "Cd9", "Itgam", "Itgax", "Lyz2", "Lirb4a", "Lgals3", "Plau", "St14", "Mmp12", "Ccl6", "Cd68", "Cxcl16", "Lgals3bp", "Gpr84", "Tspo", "Apod", "Cd86", "Trem2", "C4b", "Cd74", "Cd274", "Mpeg1"))
mgc <- AddModuleScore(
  mgc,
  gene_names,
  name = "Activation_Disease_association_score"
)

gene_names <- list(c("B2m", "Tap1", "H2-K1", "H2-D1", "H2-Q7", "H2-M3"))
mgc <- AddModuleScore(
  mgc,
  gene_names,
  name = "MHC_class_I_score"
)

gene_names <- list(c("Stat1", "Sp100", "Ifi204", "Zbp1", "Ifi44", "Isg15", "Oas2", "Oasl2", "Oasl1", "Usp18", "Irf7", "Nlrc5", "Xaf1", "Rnf213", "Ifi27l2a", "Rsad2", "Rtp4", "Ifit2", "Ifit3", "Ifit3b", "Ifit1"))
mgc <- AddModuleScore(
  mgc,
  gene_names,
  name = "Interferon_related_score"
)

gene_names <- list(c("Stat1", "Sp100", "Ifi204", "Zbp1", "Ifi44", "Isg15", "Oas2", "Oasl2", "Oasl1", "Usp18", "Irf7", "Nlrc5", "Xaf1", "Rnf213", "Ifi27l2a", "Rsad2", "Rtp4", "Ifit2", "Ifit3", "Ifit3b", "Ifit1", "B2m", "Tap1", "H2-K1", "H2-D1", "H2-Q7", "H2-M3", "Slamf9", "Cd84", "Ly9", "Ctsz", "Cst7", "Cybb", "Cd52", "Cd72", "Gpnmb", "Clec7a", "Cd9", "Itgam", "Itgax", "Lyz2", "Lirb4a", "Lgals3", "Plau", "St14", "Mmp12", "Ccl6", "Cd68", "Cxcl16", "Lgals3bp", "Gpr84", "Tspo", "Apod", "Cd86", "Trem2", "C4b", "Cd74", "Cd274", "Mpeg1"))
mgc <- AddModuleScore(
  mgc,
  gene_names,
  name = "Immunoreactive_score"
)
