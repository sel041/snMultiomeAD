set.seed(1234)
library(Seurat)
library(GenomeInfoDb)
library(GenomeInfoDbData) 
library(GenomicRanges)
library(Biobase)
library(readxl)
library(org.Mm.eg.db)
library(ggplot2) 
library(ggbeeswarm) 
library(plyr)
library(dplyr)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer) 
library(colorspace) 
library(scales) 
library(viridis)
library(VennDiagram)
options(Seurat.object.assay.version = "v5")
setwd("~/ps-renlab2/mouseAD/processed_data/gex/final")

indir <- '~/5XFAD_mouse/data/'
outdir <- '~/mouseAD/'
ct <- 'Microglia'

# temporal changes
# genotype contrasts
genotype_contrast <- 'DGE_contrasts_genotype_age'
genotype_files <- c('ADvsWT.18M.tsv', 'ADvsWT.9M.tsv', 'ADvsWT.3M.tsv')

ADvsWT.18M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.18M.tsv"))
ADvsWT.9M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.9M.tsv"))
ADvsWT.3M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.3M.tsv"))

ADvsWT.18M <- ADvsWT.18M[ADvsWT.18M$FDR < 0.05, ]
ADvsWT.9M <- ADvsWT.9M[ADvsWT.9M$FDR < 0.05, ]
ADvsWT.3M <- ADvsWT.3M[ADvsWT.3M$FDR < 0.05, ]

# up genes
ADvsWT.18M.up <- rownames(ADvsWT.18M[ADvsWT.18M$logFC > 0, ])
ADvsWT.9M.up <- rownames(ADvsWT.9M[ADvsWT.9M$logFC > 0, ])
ADvsWT.3M.up <- rownames(ADvsWT.3M[ADvsWT.3M$logFC > 0, ])

ADvsWT.early.up <- Reduce(intersect, list(ADvsWT.3M.up, ADvsWT.9M.up, ADvsWT.18M.up))
ADvsWT.mid.up <- setdiff(intersect(ADvsWT.9M.up,ADvsWT.18M.up), ADvsWT.early.up)
ADvsWT.late.up <- setdiff(ADvsWT.18M.up, c(ADvsWT.early.up, ADvsWT.mid.up))

# down genes
ADvsWT.18M.down <- rownames(ADvsWT.18M[ADvsWT.18M$logFC < 0, ])
ADvsWT.9M.down <- rownames(ADvsWT.9M[ADvsWT.9M$logFC < 0, ])
ADvsWT.3M.down <- rownames(ADvsWT.3M[ADvsWT.3M$logFC < 0, ])

ADvsWT.early.down <- Reduce(intersect, list(ADvsWT.3M.down, ADvsWT.9M.down, ADvsWT.18M.down))
ADvsWT.mid.down <- setdiff(intersect(ADvsWT.9M.down,ADvsWT.18M.down), ADvsWT.early.down)
ADvsWT.late.down <- setdiff(ADvsWT.18M.down, c(ADvsWT.early.down, ADvsWT.mid.down))

# AD- and age-associated DEGs
genotype_contrast <- 'DGE_contrasts_genotype_age'
genotype_files <- c('ADvsWT.18M.tsv', 'ADvsWT.9M.tsv', 'ADvsWT.3M.tsv')
ADvsWT.18M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.18M.tsv"))
ADvsWT.18M <- ADvsWT.18M[ADvsWT.18M$FDR < 0.05, ]

age_files <- c(#'AD.18Mvs3M.tsv', 'AD.18Mvs9M.tsv', 'AD.9Mvs3M.tsv',
               'WT.18Mvs3M.tsv', 'WT.18Mvs9M.tsv', 'WT.9Mvs3M.tsv')
WT.18Mvs3M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","WT.18Mvs3M.tsv"))
WT.18Mvs3M <- WT.18Mvs3M[WT.18Mvs3M$FDR < 0.05, ]

# up genes
ADvsWT.18M.up <- rownames(ADvsWT.18M[ADvsWT.18M$logFC > 0, ])
WT.18Mvs3M.up <- rownames(WT.18Mvs3M[WT.18Mvs3M$logFC > 0, ])
AD_Age_up <- Reduce(intersect, list(ADvsWT.18M.up, WT.18Mvs3M.up))
AD_up <- setdiff(ADvsWT.18M.up,AD_Age_up)
Age_up <- setdiff(WT.18Mvs3M.up, AD_Age_up)

# down genes
ADvsWT.18M.down <- rownames(ADvsWT.18M[ADvsWT.18M$logFC < 0, ])
WT.18Mvs3M.down <- rownames(WT.18Mvs3M[WT.18Mvs3M$logFC < 0, ])

AD_Age_down <- Reduce(intersect, list(ADvsWT.18M.down, WT.18Mvs3M.down))
AD_down <- setdiff(ADvsWT.18M.down,AD_Age_down)
Age_down <- setdiff(WT.18Mvs3M.down, AD_Age_down)
