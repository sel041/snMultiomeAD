set.seed(1234)
setwd("/tscc/projects/ps-renlab2/sel041/mouseAD/processed_data/atac/")
library(Signac)
library(Seurat)
library(JASPAR2020)

# library(Rsamtools)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(GenomicRanges)
library(Biobase)
library(readxl)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

# ggplot
library(ggplot2) 
library(ggbeeswarm) 
library(dplyr)
library(clusterProfiler)
#library(pheatmap)
#library(scCustomize)
library(pheatmap)

# color
library(RColorBrewer) 
library(colorspace) 
library(scales) 
library(viridis)

mgc <- readRDS("signac_mgc_chromvar.rds")
mgc

# Equal # of cluster
Idents(mgc) <- "final_cluster"
table(Idents(mgc))
mgc_sub <- subset(x = mgc, downsample = 619) # # of cells in the smallest cluster
mgc_sub

motif_name <- read.table("motif.names.txt", header=T)

chromvar = data.frame(mgc_sub@assays$chromvar@data)

nrow(chromvar)
length(unique(rownames(chromvar)))

all(rownames(chromvar) %in% motif_name$motif)
all(rownames(chromvar) == motif_name$motif)

motif_name_ordered <- motif_name %>%
  slice(match(rownames(chromvar), motif))

all(rownames(chromvar) == motif_name_ordered$motif)

rownames(chromvar) <- motif_name_ordered$motif.name
chromvar_filt = data.frame(chromvar)

# calculating variance
library(matrixStats)
chromvar_filt$variance = rowVars(as.matrix(chromvar_filt))
options(repr.plot.width=7, repr.plot.height=6)
hist(chromvar_filt$variance, breaks=20)

# top 10% most variable TF motifs
cutoff <- quantile(chromvar_filt$variance, 0.9, na.rm = TRUE)
cutoff

# filtering for top variable motifs
Top_variance = subset(chromvar_filt, chromvar_filt$variance >= cutoff)
nrow(Top_variance)

#subsetting chromvar data frame for top varaible motifs
chromvar_top = subset(chromvar_filt, rownames(chromvar_filt) %in% rownames(Top_variance))
#removing column having variance values from subseted chromvar data frame
chromvar_top = chromvar_top[,-23666]
#calculating motif correlation
chromvar_cor=cor(t(chromvar_top))
#hirarchical clustering with pheatmap
options(repr.plot.width=10, repr.plot.height=10)
p <- pheatmap(chromvar_cor,color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                "RdBu")))(100), border_color = NA, fontsize = 6)
p

# plot a violin plot 
options(repr.plot.width=3.2, repr.plot.height=4)
motif_to_plot = "IRF7"
vln_plot <- VlnPlot(mgc_sub, features = motif_name[motif_name$motif.name == motif_to_plot, ]$motif, 
                    cols = c("#777ED8", "#d486d1", "#45ddce"), 
                    pt.size = 0)+NoLegend() # Disable the default dot plot to avoid clutter

vln_plot <- vln_plot +
  geom_boxplot(width = 0.2, 
               color = "black", 
               alpha = 0.2, 
               outlier.shape = NA) + # Remove outliers from the box plot
  stat_summary(fun = mean, 
               geom = "line", 
               shape = 18, 
               size = 3, 
               color = "black") + # Add mean points
  theme_classic() +
  theme(
    text = element_text(color = "black", size=12),       
    axis.title = element_text(color = "black", size=12), 
    axis.text = element_text(color = "black", size=12), 
    plot.title = element_text(color = "black", size=12) 
  ) +
  labs(title = c(),
       x = "Group",
       y = paste0(motif_to_plot, " Chromvar Activity"))
vln_plot
