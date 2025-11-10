## pdac subtype
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(data.table)
library(tidyverse)
library(dplyr)


# load object
load("./fitteds_public_2019-02-12.Rdata")
source("purIST_function.R")

Mx <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  row.names = "Geneid",
  check.names = FALSE
)

Gene_length <- Mx %>% 
  #select(Geneid, Length) %>% 
  #remove_rownames() %>% 
  #column_to_rownames(var = "Geneid") %>% 
  pull(Length)
CountTable <- Mx %>% 
  select(-c(2:6)) %>% 
  #remove_rownames() %>% 
  #column_to_rownames(var = "Geneid") %>% 
  rename_with(~ str_replace_all(., "-.*", ""))

biomarkers = c("GPR87", "REG4", 
               "KRT6A", "ANXA10",
               "BCAR3", "GATA6",
               "PTGES", "CLDN18",
               "ITGA3", "LGALS4",
               "C16orf74", "DDC", 
               "S100A2", "SLC40A1", 
               "KRT5", "CLRN3")
  
### TPM Normalization
# Step 1: Calculate Reads Per Kilobase (RPK)
rpk <- CountTable / (Gene_length / 1000)

# Step 2: Calculate per-sample scaling factor (sum of RPKs)
scaling_factors <- colSums(rpk)

# Step 3: Calculate TPM
CountTable_TPM <- t(t(rpk) / scaling_factors * 1e6)

CountTable_TPM_mini <- CountTable_TPM %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "row_names") %>% 
  filter(row_names %in% biomarkers) %>% 
  column_to_rownames(var = "row_names")

pred <- apply_classifier(data = CountTable_TPM_mini, classifier = classifs[[1]])

fwrite(pred, file = snakemake@output[["out_pred"]], quote = F, sep = "\t", row.names = T, col.names = T)


### heatmap
df_colData <- data.frame(Subtype = pred$Subtype_graded)
rownames(df_colData) <- rownames(pred)

dds <- DESeqDataSetFromMatrix(countData = CountTable[biomarkers,], colData = df_colData, design = ~1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
vds <- DESeq2::varianceStabilizingTransformation(assay(dds))
img <- vds %>% 
  pheatmap::pheatmap(show_colnames = T, show_rownames = T,
                     scale = "row",
                     #cutree_cols = 8,
                     clustering_method = "complete",
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     annotation_col = df_colData,
                     #annotation_colors = my_Anno_cols,
                     cluster_cols = T, 
                     angle_col = 315)

ggsave(snakemake@output[["out_heatmap"]], img)
#png("./img/p1713_pheatmap_pdac_subtype_2.png")
#img
#dev.off()