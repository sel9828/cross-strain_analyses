# Project: Microalgae Cross-strain Growth Medium Reuse Experiment
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# This code carries out differential abundance analyses on bacteria OTU datasets.
# Notes on Deseq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Input data and metadata are available on Figshare: https://doi.org/10.6084/m9.figshare.7831913.v2

# Load packages
library(DESeq2) # version 1.14.1
library(tidyverse)
library(gplots)

# Requirements for DESeq dataset:
# 1. cts = Counts matrix in absolute numbers (integers, not normalized) (OTUs as rows, samples as columns)
# 2. coldata = Table of column data; first column row names must match columns names of cts; additional columns describe the sample (e.g., treatment, replicate)
# 3. design = a column of coldata that is being used to compare the samples (e.g., treatment)

# Counts matrix (including Spirochaete OTUs)
otutable <- read.delim("input_data/Data3_Absolute_OTUs.txt", skip = 1, header = T, row.names = 1)  # OTU table

otutable_no_inoc <- otutable[ , -c(1, 8, 15)] # remove inoculum samples, since only comparing treatment and control on final day

cts <- round(otutable_no_inoc[ , !(colnames(otutable_no_inoc) %in% c("taxonomy"))]) # remove taxonomy column, leaving only sample names as columns; round numbers to whole integers

# Split cts matrix by algae 
cts_C323 <- cts[,1:6]   # Staurosira sp. C323
cts_D046 <- cts[,7:12]  # Chlorella sp. D046
cts_Navi <- cts[,13:18] # Navicula sp. SFP

# Remove OTUs that have 0 counts in all samples
cts_C323 <- cts_C323[rowSums(cts_C323)!=0, ]
cts_D046 <- cts_D046[rowSums(cts_D046)!=0, ]
cts_Navi <- cts_Navi[rowSums(cts_Navi)!=0, ]

# Column data 
reps <- c("A", "B", "C", "D", "E", "F") # biological replicates
treatment <- c(rep("Fresh", 3), rep("Reused", 3))

coldata_C323 <- data.frame(reps, treatment)
coldata_D046 <- data.frame(reps, treatment)
coldata_Navi <- data.frame(reps, treatment)

rownames(coldata_C323) <- colnames(cts_C323)
rownames(coldata_D046) <- colnames(cts_D046)
rownames(coldata_Navi) <- colnames(cts_Navi)

# DESeq dataset
dds_C323 <- DESeqDataSetFromMatrix(countData = cts_C323, 
                              colData = coldata_C323,
                              design = ~ treatment)

dds_D046 <- DESeqDataSetFromMatrix(countData = cts_D046, 
                                  colData = coldata_D046,
                                  design = ~ treatment)

dds_Navi <- DESeqDataSetFromMatrix(countData = cts_Navi, 
                                  colData = coldata_Navi,
                                  design = ~ treatment)

# Specify the control treatment
dds_C323$treatment <- relevel(dds_C323$treatment, ref = "Fresh")
dds_D046$treatment <- relevel(dds_D046$treatment, ref = "Fresh")
dds_Navi$treatment <- relevel(dds_Navi$treatment, ref = "Fresh")

# Differential analysis. 
  # MAP refers to methods that uses LC shrinkage, and calculates a p-value for the shrunken LFC.
   # MAP does Log Fold Change shrinkage (*note: Need DESeq2 version 1.16 for lfcShrink function, requires R version 3.4)
      # From the beginner's guide: "When count values are too low to allow an accurate estimate of the LFC, 
      # the value is “shrunken” towards zero to avoid that these values, which otherwise would frequently be 
      # unrealistically large, dominate the top-ranked log fold changes."

  # Decided to use MLE for p-values/adjp and log fold change (since shrinkage MAP may not be relevant to this type of data)

# MLE
DESeq_C323_MLE <- DESeq(dds_C323, modelMatrixType="standard", betaPrior=FALSE)
DESeq_D046_MLE <- DESeq(dds_D046, modelMatrixType="standard", betaPrior=FALSE)
DESeq_Navi_MLE <- DESeq(dds_Navi, modelMatrixType="standard", betaPrior=FALSE)

# MAP
# DESeq_C323_MAP <- DESeq(dds_C323, modelMatrixType="standard", betaPrior=TRUE)
# DESeq_D046_MAP <- DESeq(dds_D046, modelMatrixType="standard", betaPrior=TRUE)
# DESeq_Navi_MAP <- DESeq(dds_Navi, modelMatrixType="standard", betaPrior=TRUE)

# Differential analysis results
  # Change the alpha value (false discovery rate cutoff for p-adjusted) for independent filtering step 
    # to 0.05 since this is what I'll be using for determining significance
  # padj = adjusted p value, is the fraction of false positives if one called significant the value of this OTU's p-value or lower

a <- 0.05

# MLE results 
res_C323_MLE <- results(DESeq_C323_MLE, name="treatment_Reused_vs_Fresh", alpha = a)
res_D046_MLE <- results(DESeq_D046_MLE, name="treatment_Reused_vs_Fresh", alpha = a)
res_Navi_MLE <- results(DESeq_Navi_MLE, name="treatment_Reused_vs_Fresh", alpha = a)

# MAP results 
# res_C323_MAP <- results(DESeq_C323_MAP, name="treatment_Recycled_vs_Fresh", alpha = a)
# res_D046_MAP <- results(DESeq_D046_MAP, name="treatment_Recycled_vs_Fresh", alpha = a)
# res_Navi_MAP <- results(DESeq_Navi_MAP, name="treatment_Recycled_vs_Fresh", alpha = a)

# Summaries of adjp < a
summary(res_C323_MLE)
summary(res_D046_MLE)
summary(res_Navi_MLE)

# Export results 
write.csv(as.data.frame(res_C323_MLE), file="DESeq2_results/C323_DESeq2_results.csv")
write.csv(as.data.frame(res_D046_MLE), file="DESeq2_results/D046_DESeq2_results.csv")
write.csv(as.data.frame(res_Navi_MLE), file="DESeq2_results/Navi_DESeq2_results.csv")

# Tables of sig results 
res_C323_sig <- as.data.frame(subset(res_C323_MLE, padj < a))
res_D046_sig <- subset(res_D046_MLE, padj < a)
res_Navi_sig <- subset(res_Navi_MLE, padj < a)


# # Merge the log2FoldChange results from the three algae to make a heatmap (not published)

# res_C323_log2 <- as.data.frame(res_C323_MLE)
# res_C323_log2$OTU <- rownames(res_C323_log2)
# res_C323_log2 %>%
#   select(OTU, C323_log2 = log2FoldChange) -> C323_log2
# 
# res_D046_log2 <- as.data.frame(res_D046_MLE)
# res_D046_log2$OTU <- rownames(res_D046_log2)
# res_D046_log2 %>%
#   select(OTU, D046_log2 = log2FoldChange) -> D046_log2
# 
# res_Navi_log2 <- as.data.frame(res_Navi_MLE)
# res_Navi_log2$OTU <- rownames(res_Navi_log2)
# res_Navi_log2 %>%
#   select(OTU, Navi_log2 = log2FoldChange) -> Navi_log2
# 
# D046_log2 %>%
#   full_join(Navi_log2, by = "OTU") %>%
#   full_join(C323_log2, by = "OTU") -> all_algae_log2
# 
# rownames(all_algae_log2) <- all_algae_log2$OTU
# 
# library(RColorBrewer)
# my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(100)
# 
# png("figures/deseq2_log2_heatmap.png",
#     width = 7*300,
#     height = 6.5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)
# heatmap.2(as.matrix(all_algae_log2[  ,2:4]),
#           density.info = "none", trace = "none", dendrogram = "none",
#           na.color = "white", Colv = FALSE, Rowv = FALSE, col = my_palette,
#           key = T, labCol = NA, key.xlab = "Log2 Fold Change", key.title = "", keysize = 1)
# dev.off()


