# Marker Overlap Analysis
# Code adapted by Caroline He from https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/basic_tutorial.Rmd
#2023

## 1. load packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)

## 2. Set random seed for reproducibility
set.seed(12345)

## 3. Set working directory
setwd("~/new_frontal_official/hdWGCNA/Ast-subset")

## 4. Load the seurat object that is made at the end of the hdWGCNA-module_identification.R script
seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/Ast-subset/hdWGCNA-seurat_obj-basic-end.rds")

## 5. Load the list of marker genes associated with each astrocyte cluster
Ast_clusters_pos <- read_csv("Ast_clusters_pos.csv")
Ast_clusters_neg <- read_csv("Ast_clusters_neg.csv")

markers <- as.data.frame(Ast_clusters_pos)
  
## 6. compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 0.25 # log fold change cutoff for overlap analysis
)

write.csv(overlap_df, "Ast-cluster_overlap_UP.csv")

## 7. Visualize the overlap and associated statistics with the gene list comparison
plot_list <- OverlapBarPlot(overlap_df)

wrap_plots(plot_list, ncol=3)

OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & positive cluster markers')



