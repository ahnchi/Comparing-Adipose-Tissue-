#Required packages
library(stringr)
library(PhosR)
library(scales)
library(devEMF)
library(limma)
library(stats)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(viridis)
library(tidyverse)
library(ComplexHeatmap)
library(dendextend)
library(devtools)
library(edgeR)
library(grDevices)

#Set working directory
setwd("")
#Import data
#"phospho_imputed.csv" is a phosphoproteomics data that is 1) normalized by median level for each phosphopeptide, 
# and 2) imputed by using 'impseq' method from NAguideR. http://www.omicsolution.org/wukong/NAguideR/
peptide_imputed <- read.csv("phospho_imputed.csv", row.names = 1)

# Collapse phosphosites
id1 <- sapply(strsplit(rownames(peptide_imputed), "~"), function(x){paste(toupper(x[2]), x[3], x[4], sep=";")})
phospho.site1 <- phosCollapse(peptide_imputed, id = id1, stat = apply(abs(peptide_imputed), 1, max, na.rm=TRUE), by = "max")

##phospho.site1 is median normalized, imputed by NAguideR. Let's call this as phospho.filtered1
phospho.filtered1 <- phospho.site1

num_rows_with_negative <- sum(rowSums(phospho.filtered1 < 0) >= 1)
print(num_rows_with_negative)
# Remove rows with at least one negative value
phospho.filtered1 <- phospho.filtered1[rowSums(phospho.filtered1 < 0) == 0, ]

#edgeR--------------------------------------------------------------------------
group <- factor(rep(c("EX", "SED"), each = 8))
# Create a pairing factor for the paired samples
pairing <- factor(rep(1:8, times = 2))
# Use relevel to set 'SED' as the reference level
group_relevel <- relevel(group, ref = "SED")
dge <- DGEList(counts = mat, group = group_relevel)
# 3. Estimate dispersion
dge <- estimateDisp(dge, design = model.matrix(~group+pairing))
plotBCV(dge)

# Perform differential expression analysis
et <- exactTest(dge)
tags <- topTags(et, n = Inf)

# extract the results table and add gene information
results_phospho <- tags$table

#Pathway analysis (Reactome)----------------------------------------------------
Tc.gene <- phosCollapse(results_phospho, id=gsub(";.+", "", rownames(results_phospho)), 
                        stat=apply(abs(results_phospho), 1, max), by = "max")
Tc.gene <- as.matrix(Tc.gene[, c(1)])

#Load required packages
library(org.Hs.eg.db)
library(reactome.db)
library(annotate)
#Prep Reactome db.
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]
pathways = pathways[which(grepl("Homo sapiens: ", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

data(PM)
data(Pathways)

#Rank based pathway analysis 
path_greater <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "greater") #Enriched in EX.
path_less <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                   annotation=pathways, 
                                   alter = "less") #Enriched in SED.
