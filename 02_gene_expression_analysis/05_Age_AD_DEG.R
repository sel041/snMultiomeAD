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

# AD temporal changes 
# genotype contrasts
genotype_contrast <- 'DGE_contrasts_genotype_age'
genotype_files <- c('ADvsWT.18M.tsv', 'ADvsWT.9M.tsv', 'ADvsWT.3M.tsv')

ADvsWT.18M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.18M.tsv"))
ADvsWT.9M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.9M.tsv"))
ADvsWT.3M <- read.table(paste0(indir,genotype_contrast,"/",ct,"/","ADvsWT.3M.tsv"))

ADvsWT.18M <- ADvsWT.18M[ADvsWT.18M$FDR < 0.05, ]
ADvsWT.9M <- ADvsWT.9M[ADvsWT.9M$FDR < 0.05, ]
ADvsWT.3M <- ADvsWT.3M[ADvsWT.3M$FDR < 0.05, ]

ADvsWT.18M <- ADvsWT.18M[!grepl("^Gm|^Rpl|^Rps|^Mrps|^Mrpl", rownames(ADvsWT.18M)), ] 
ADvsWT.9M <- ADvsWT.9M[!grepl("^Gm|^Rpl|^Rps|^Mrps|^Mrpl", rownames(ADvsWT.9M)), ]
ADvsWT.3M <- ADvsWT.3M[!grepl("^Gm|^Rpl|^Rps|^Mrps|^Mrpl", rownames(ADvsWT.3M)), ]
