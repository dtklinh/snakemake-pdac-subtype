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
  remove_rownames() %>% 
  column_to_rownames(var = "Geneid") %>% 
  rename_with(~ str_replace_all(., "-.*", ""))