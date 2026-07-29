# soloTE heatmap plot
conditions <- c("WT", "AD")
timepoints <- c("3M", "9M", "18M")
for (cond in conditions){
  for (tp in timepoints){
    group_name <- paste0(cond,'_',tp)
    path2matrix = paste0('SoloTE_5xFAD/Microglia_',cond,'_',tp,'_SoloTE_output/Microglia_',cond,'_',tp,'_familytes_MATRIX/')
    solote_matrix <- ReadMtx(paste0(path2matrix,"matrix.mtx"),paste0(path2matrix,"barcodes.tsv"),paste0(path2matrix,"features.tsv"),feature.column=1)
    obj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
    obj$condition <- cond 
    obj$timepoint <- tp
    obj$group <- group_name  
    assign (group_name, obj)
  }
}

# Merge all objects
seurat_merged <- merge(
  x    = WT_3M,
  y    = list(WT_9M, WT_18M, AD_3M, AD_9M, AD_18M),
  add.cell.ids = c("WT_3M", "WT_9M", "WT_18M", "AD_3M", "AD_9M", "AD_18M"),
  merge.data   = TRUE
)

# Verify
table(seurat_merged$group)
table(seurat_merged$condition, seurat_merged$timepoint)

# Normalize
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, features = all.genes)

# Join the layers 
seurat_merged <- JoinLayers(seurat_merged)
Idents(seurat_merged) <- "group"

# Fold change compared to 3-month-old WT
wt9m <- FindMarkers(
  seurat_merged,
  ident.1   = "WT_9M",    
  ident.2   = "WT_3M",
  only.pos  = FALSE,      
  min.pct   = 0,        
  logfc.threshold = 0  
)
wt9m$gene <- rownames(wt9m)

wt18m <- FindMarkers(
  seurat_merged,
  ident.1   = "WT_18M",    
  ident.2   = "WT_3M",
  only.pos  = FALSE,      
  min.pct   = 0,        
  logfc.threshold = 0  
)
wt18m$gene <- rownames(wt18m)

ad3m <- FindMarkers(
  seurat_merged,
  ident.1   = "AD_3M",    
  ident.2   = "WT_3M",
  only.pos  = FALSE,      
  min.pct   = 0,        
  logfc.threshold = 0  
)
ad3m$gene <- rownames(ad3m)

ad9m <- FindMarkers(
  seurat_merged,
  ident.1   = "AD_9M",    
  ident.2   = "WT_3M",
  only.pos  = FALSE,      
  min.pct   = 0,        
  logfc.threshold = 0  
)
ad9m$gene <- rownames(ad9m)

ad18m <- FindMarkers(
  seurat_merged,
  ident.1   = "AD_18M",    
  ident.2   = "WT_3M",
  only.pos  = FALSE,      
  min.pct   = 0,        
  logfc.threshold = 0  
)
wt18m$gene <- rownames(wt18m)

wt3m <- wt9m
wt3m$pct.1 <- 1
wt3m$pct.2 <- 1
wt3m$p_val <- 1
wt3m$p_val_adj <- 1
wt3m$avg_log2FC <- 0
head(wt3m)

wt3m$gene <- gsub("SoloTE-", "", wt3m$gene)
wt9m$gene <- gsub("SoloTE-", "", wt9m$gene)
wt18m$gene <- gsub("SoloTE-", "", wt18m$gene)

ad3m$gene <- gsub("SoloTE-", "", ad3m$gene)
ad9m$gene <- gsub("SoloTE-", "", ad9m$gene)
ad18m$gene <- gsub("SoloTE-", "", ad18m$gene)

library(dplyr)
library(ggplot2)
wt3m$group  <- "wt3m"
wt9m$group  <- "wt9m"
wt18m$group <- "wt18m"
ad3m$group  <- "ad3m"
ad9m$group  <- "ad9m"
ad18m$group <- "ad18m"

df <- bind_rows(wt3m, wt9m, wt18m, ad3m, ad9m, ad18m)
genes_to_remove <- unique(df$gene[df$pct.2 == 0])
df <- df %>%
  filter(!gene %in% genes_to_remove)
df <- df %>%
  mutate(sig = case_when(
    p_val_adj < 0.0001 ~ "***",
    p_val_adj < 0.001  ~ "**",
    p_val_adj < 0.01   ~ "*",
    TRUE ~ ""
  ))
df$group <- factor(df$group,
                   levels = c("wt3m","wt9m","wt18m",
                              "ad3m","ad9m","ad18m"))
df$gene <- factor(df$gene, levels = unique(df$gene))
gene_levels <- levels(df$gene)

# Plot heatmap
p <- ggplot(df, aes(x = gene, y = group, fill = avg_log2FC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sig), size = 2.5) +
  
  scale_fill_distiller(
    palette = "RdYlBu",
    direction = -1,
    limits = c(-4, 5)   
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Gene", y = NULL, fill = "avg_log2FC")
p
