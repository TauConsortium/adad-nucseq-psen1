#https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/basic_tutorial.Rmd

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)

# set random seed for reproducibility
set.seed(12345)

setwd("~/new_frontal_official/hdWGCNA/Ast-subset")

seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/Ast-subset/hdWGCNA-seurat_obj-basic-end.rds")

# compute cell-type marker genes with Seurat:
Ast_clusters_pos <- read_csv("Ast_clusters_pos.csv")
Ast_clusters_neg <- read_csv("Ast_clusters_neg.csv")

markers <- as.data.frame(Ast_clusters_pos)
  
# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 0.25 # log fold change cutoff for overlap analysis
)

#setwd("~/new_frontal_official/hdWGCNA/Ast-subset/")
#write.csv(overlap_df, "Ast-cluster_overlap_DOWN.csv")

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=3)

OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & positive cluster markers')



