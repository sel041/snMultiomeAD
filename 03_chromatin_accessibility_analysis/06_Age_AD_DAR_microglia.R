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
setwd("~/mouseAD")

indir <- '~/5XFAD_mouse/data/'
outdir <- '~/mouseAD/'
ct <- 'Microglia'

# temporal changes
# genotype contrasts
genotype_contrast <- 'DAR_contrasts_genotype_age'
genotype_files <- c('ADvsWT.18M.tsv', 'ADvsWT.9M.tsv', 'ADvsWT.3M.tsv')

ADvsWT.18M.DAR <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.18M.tsv"))
ADvsWT.9M.DAR <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.9M.tsv"))
ADvsWT.3M.DAR <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.3M.tsv"))

ADvsWT.18M.DAR <- ADvsWT.18M.DAR[ADvsWT.18M.DAR$FDR < 0.05, ]
ADvsWT.9M.DAR <- ADvsWT.9M.DAR[ADvsWT.9M.DAR$FDR < 0.05, ]
ADvsWT.3M.DAR <- ADvsWT.3M.DAR[ADvsWT.3M.DAR$FDR < 0.05, ]

# up
ADvsWT.18M.DAR.up <- rownames(ADvsWT.18M.DAR[ADvsWT.18M.DAR$logFC > 0, ])
ADvsWT.9M.DAR.up <- rownames(ADvsWT.9M.DAR[ADvsWT.9M.DAR$logFC > 0, ])
ADvsWT.3M.DAR.up <- rownames(ADvsWT.3M.DAR[ADvsWT.3M.DAR$logFC > 0, ])

# down
ADvsWT.18M.DAR.down <- rownames(ADvsWT.18M.DAR[ADvsWT.18M.DAR$logFC < 0, ])
ADvsWT.9M.DAR.down <- rownames(ADvsWT.9M.DAR[ADvsWT.9M.DAR$logFC < 0, ])
ADvsWT.3M.DAR.down <- rownames(ADvsWT.3M.DAR[ADvsWT.3M.DAR$logFC < 0, ])

# age & genotype contrasts
genotype_contrast <- 'DAR_contrasts_genotype_age'
genotype_files <- c('ADvsWT.18M.tsv', 'ADvsWT.9M.tsv', 'ADvsWT.3M.tsv')

ADvsWT.18M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.18M.tsv"))
ADvsWT.18M <- ADvsWT.18M[ADvsWT.18M$FDR < 0.05, ]

age_contrast <- 'DAR_contrasts_genotype_age'
age_files <- c(#'AD.18Mvs3M.tsv', 'AD.18Mvs9M.tsv', 'AD.9Mvs3M.tsv',
               'WT.18Mvs3M.tsv', 'WT.18Mvs9M.tsv', 'WT.9Mvs3M.tsv')

WT.18Mvs3M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","WT.18Mvs3M.tsv"))
WT.18Mvs3M <- WT.18Mvs3M[WT.18Mvs3M$FDR < 0.05, ]

# up DARs
ADvsWT.18M.up <- rownames(ADvsWT.18M[ADvsWT.18M$logFC > 0, ])
WT.18Mvs3M.up <- rownames(WT.18Mvs3M[WT.18Mvs3M$logFC > 0, ])
AD_Age_up <- Reduce(intersect, list(ADvsWT.18M.up, WT.18Mvs3M.up))
AD_up <- setdiff(ADvsWT.18M.up,AD_Age_up)
Age_up <- setdiff(WT.18Mvs3M.up, AD_Age_up)

# down DARs
ADvsWT.18M.down <- rownames(ADvsWT.18M[ADvsWT.18M$logFC < 0, ])
WT.18Mvs3M.down <- rownames(WT.18Mvs3M[WT.18Mvs3M$logFC < 0, ])
AD_Age_down <- Reduce(intersect, list(ADvsWT.18M.down, WT.18Mvs3M.down))
AD_down <- setdiff(ADvsWT.18M.down,AD_Age_down)
Age_down <- setdiff(WT.18Mvs3M.down, AD_Age_down)




