library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# heatmaps TE family compared with WT 3M
counts <- read.table('microglia_geno_age_counts_union.tsv',header=T)
cpm <- apply(counts, 2, function(x) {(x/sum(x)*1E6)})
cpm <- as.data.frame(log2(cpm+1))
cpm$peaks <- row.names(counts)

union <- read.table('../macs2_new/q_05/final_peaks/mouse_ad_union.bed')
union <- union %>% unite("peaks", V1, V2, V3, sep = "-")
cpm <- join(cpm, union, by = 'peaks')

family <- read.table('mouse_ad_union_peaks_te_family.tsv',header=F, col.names=c('V4','TE_Family'))
cpm <- join(cpm, family, by = 'V4')

cpm <- cpm[!grepl("\\?", cpm$TE_Family), ]
cpm <- na.omit(cpm)

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$WT_9M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_9M_WT <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$WT_18M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_18M_WT <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$X5XFAD_3M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_3M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$X5XFAD_9M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_9M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$X5XFAD_18M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_18M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

p_18M_WT$TE_Family <- NULL
p_3M_5xfad$TE_Family <- NULL
p_9M_5xfad$TE_Family <- NULL
p_18M_5xfad$TE_Family <- NULL

pvalues <- cbind(as.numeric('1'), p_9M_WT, p_18M_WT, p_3M_5xfad, p_9M_5xfad, p_18M_5xfad)
colnames(pvalues) <- c('p_WT_3M', 'TE_Family', 'p_WT_9M', 'p_WT_18M','p_5XFAD_3M','p_5XFAD_9M','p_5XFAD_18M')
cpm$WT_9M <- log2((cpm$WT_9M+0.1)/(cpm$WT_3M+0.1))
cpm$WT_18M <- log2((cpm$WT_18M+0.1)/(cpm$WT_3M+0.1))
cpm$X5XFAD_3M <- log2((cpm$X5XFAD_3M+0.1)/(cpm$WT_3M+0.1))
cpm$X5XFAD_9M <- log2((cpm$X5XFAD_9M+0.1)/(cpm$WT_3M+0.1))
cpm$X5XFAD_18M <- log2((cpm$X5XFAD_18M+0.1)/(cpm$WT_3M+0.1))
cpm$WT_3M <- 0

cpm$peaks <- NULL
cpm$V4 <- NULL
cpm <- cpm[!grepl("\\?", cpm$TE_Family), ]
cpm_avg <- cpm %>%
  group_by(TE_Family) %>%
  summarise(across(1:6, mean, na.rm = TRUE))

cpm_avg <- as.data.frame(cpm_avg)

cpm_avg <- join(cpm_avg, pvalues, by = 'TE_Family')
rownames(cpm_avg) <- cpm_avg$TE_Family
cpm_avg$TE_Family <- NULL

anno <- matrix("", nrow = nrow(cpm_avg[,1:6]), ncol = ncol(cpm_avg[,1:6]))
colnames(anno) <- colnames(cpm_avg[,1:6])
rownames(anno) <- rownames(cpm_avg[,1:6])
anno[cpm_avg[,7:12]<0.01] <- "*"
anno[cpm_avg[,7:12]<0.001] <- "**"
anno[cpm_avg[,7:12]<0.0001] <- "***"

options(repr.plot.width = 4.2, repr.plot.height = 12, repr.plot.res = 200)
hm <- pheatmap(mat=cpm_avg[,1:6], 
         scale = 'none', 
         border_color = 'black', 
         cluster_cols = F, 
         cluster_rows = T, 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         fontsize_row = 12,
         fontsize_col = 12,
         display_numbers = anno,
         fontsize_number = 10,
         main = "         Chrom Access at TE Families")

hm


# heatmaps TE family 5xFAD / WT
counts <- read.table('microglia_geno_age_counts_union.tsv',header=T)
cpm <- apply(counts, 2, function(x) {(x/sum(x)*1E6)})
cpm <- as.data.frame(log2(cpm+1))
cpm$peaks <- row.names(counts)

union <- read.table('../macs2_new/q_05/final_peaks/mouse_ad_union.bed')
union <- union %>% unite("peaks", V1, V2, V3, sep = "-")
cpm <- join(cpm, union, by = 'peaks')

family <- read.table('mouse_ad_union_peaks_te_family.tsv',header=F, col.names=c('V4','TE_Family'))
cpm <- join(cpm, family, by = 'V4')
cpm <- cpm[!grepl("\\?", cpm$TE_Family), ]
cpm <- na.omit(cpm)

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_3M, cpm_subset$X5XFAD_3M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_3M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_9M, cpm_subset$X5XFAD_9M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_9M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

# Initialize a list to store the p-values for each TE_Family
p_values <- list()

# Loop through each unique TE_Family and perform the Wilcoxon test
for (te_family in unique(cpm$TE_Family)) {
  # Filter the data for the current TE_Family
  cpm_subset <- cpm %>% filter(TE_Family == te_family)
  
  # Perform the Wilcoxon test between WT_3M and WT_9M
  if (nrow(cpm_subset) > 1) {  # Ensure there are enough rows to perform the test
    test_result <- wilcox.test(cpm_subset$WT_18M, cpm_subset$X5XFAD_18M, paired = TRUE)
    # Store the p-value in the list with the TE_Family as the name
    p_values[[te_family]] <- test_result$p.value
  } else {
    # If there's not enough data, store NA for this TE_Family
    p_values[[te_family]] <- NA
  }
}

# Convert the list of p-values to a data frame for easier viewing
p_18M_5xfad <- data.frame(TE_Family = names(p_values), p_value = unlist(p_values))

p_3M_5xfad$TE_Family <- NULL
p_9M_5xfad$TE_Family <- NULL
p_18M_5xfad$TE_Family <- NULL

pvalues <- cbind(p_3M_5xfad, p_9M_5xfad, p_18M_5xfad)
pvalues$TE_Family <- rownames(pvalues)
colnames(pvalues) <- c('p_5XFAD_3M','p_5XFAD_9M','p_5XFAD_18M','TE_Family')
cpm$X5XFAD_3M <- log2((cpm$X5XFAD_3M+0.1)/(cpm$WT_3M+0.1))
cpm$X5XFAD_9M <- log2((cpm$X5XFAD_9M+0.1)/(cpm$WT_9M+0.1))
cpm$X5XFAD_18M <- log2((cpm$X5XFAD_18M+0.1)/(cpm$WT_18M+0.1))

cpm[,1:3] <- NULL
cpm$peaks <- NULL
cpm$V4 <- NULL
cpm <- cpm[!grepl("\\?", cpm$TE_Family), ]
cpm_avg <- cpm %>%
  group_by(TE_Family) %>%
  summarise(across(1:3, mean, na.rm = TRUE))

cpm_avg <- as.data.frame(cpm_avg)
cpm_avg <- join(cpm_avg, pvalues, by = 'TE_Family')
rownames(cpm_avg) <- cpm_avg$TE_Family
cpm_avg$TE_Family <- NULL

anno <- matrix("", nrow = nrow(cpm_avg[,1:3]), ncol = ncol(cpm_avg[,1:3]))
colnames(anno) <- colnames(cpm_avg[,1:3])
rownames(anno) <- rownames(cpm_avg[,1:3])
anno[cpm_avg[,4:6]<0.01] <- "*"
anno[cpm_avg[,4:6]<0.001] <- "**"
anno[cpm_avg[,4:6]<0.0001] <- "***"

hm <- pheatmap(mat=cpm_avg[,1:3],
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(60),
         breaks = seq(-0.5, 0.5, 1/60), 
         scale = 'none', 
         border_color = 'black', 
         cluster_cols = F, 
         cluster_rows = T, 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         fontsize_row = 12,
         fontsize_col = 12,
         display_numbers = anno,
         fontsize_number = 10,
         main = "         Chrom Access at TE Families")

hm

