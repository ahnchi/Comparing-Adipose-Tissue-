# Used packages
library(stringr)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)

#Set working directory
setwd("")

#Read proteomics data
#SED and EX are arranged in paired order. E.g., SED_1 is matched with EX_1 ... SED_8 is matched with EX_8
df <- read.csv("proteomics.csv", row.names = 1)

#Normalization by median level--------------------------------------------------
medians <- apply(df[, 5:20], 2, median)
max_val2 <- max(medians)
# calculate the normalizing factor for each column
norm_factors2 <- max_val2/medians
norm_factors2 <- unname(norm_factors2)

###so now applying to df
df_norm_median <- sweep(df[,5:20], 2, norm_factors2, '*')  ##df_norm_median is normalized data that only contains proteomics matrix.
df_norm_m <- cbind(df[,1:4], df_norm_median)   ##df_norm_m is now ready for stats test. 

###edgeR------------------------------------------------------------------------
mat <- as.matrix(df_norm_median)
group <- factor(rep(c("SED", "EX"), each = 8))
# Create a pairing factor for the paired samples
pairing <- factor(rep(1:8, times = 2))
# Use relevel to set 'SED' as the reference level
group_relevel <- relevel(group, ref = "SED")

dge <- DGEList(counts = mat, group = group_relevel)
dge <- estimateDisp(dge, design = model.matrix(~group+pairing))
plotBCV(dge)

#Perform differential expression analysis
et <- exactTest(dge)
tags <- topTags(et, n = Inf)
results <- tags$table
matched_df_norm_m <- df_norm_m[row.names(results), c("Accession", "Description", "GeneID")]
merged_results <- merge(matched_df_norm_m, results, by = 0, all.x = TRUE)

#Clean up
merged_results <- merged_results[,-1] 
