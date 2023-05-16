# Weighted Gene Coexpression Network Analysis (WGCNA)
# Code adapted by Caroline He from https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/basic_tutorial.Rmd
# 2023

# INPUT: Describe input
load("~/new_frontal_official/brainFRONTAL.noMT.RData")


# STEP 1: PREPROCESSING

# 1. Load libraries:
## single-cell analysis package
library(Seurat)
## plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
## co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
## using the cowplot theme for ggplot
theme_set(theme_cowplot())
## set random seed for reproducibility
set.seed(12345)

# 2. Set local directory
setwd("~/new_frontal_official/hdWGCNA/Ast-subset")

# 3. Build Seurat object
# This is how the file had been processed before:
#seurat_obj<-NormalizeData(seurat_obj)
#seurat_obj<-FindVariableFeatures(seurat_obj)
#seurat_obj<-ScaleData(seurat_obj)
#seurat_obj<-RunPCA(seurat_obj, features=VariableFeatures(seurat_obj))
#seurat_obj<-RunUMAP(seurat_obj, dims=1:20)
#seurat_obj<-FindNeighbors(seurat_obj)
#seurat_obj<-FindClusters(seurat_obj)

seurat_obj<-BrainFRONTAL.noMT  #Miochondrial genes have been removed
seurat_obj<- subset(x = seurat_obj, subset = Cluster.id == "Astrocytes") #Retaining only astrocytes from the Seurat object
DefaultAssay(seurat_obj)<- 'RNA'

# 4. Load Seurat object for WGCNA 
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "BrainFP" # the name of the hdWGCNA experiment
)

# 5. Aggregate cells into 'meta-cells'. Here we group 15 single cells into a "meta-cell"
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Patient", "Cluster.id"), # specify the columns in seurat_obj@meta.data to group by
  k = 15, # nearest-neighbors parameter, use 15 cells
  ident.group = 'Cluster.id' # set the Idents of the metacell seurat object
)

# 6. Normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- seurat_obj %>%
  NormalizeMetacells() %>%
  ScaleMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunPCAMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunHarmonyMetacells(group.by.vars='Patient') %>%
  RunUMAPMetacells(reduction='harmony', dims=1:15)
p1 <- DimPlotMetacells(seurat_obj, group.by='Cluster.id') + umap_theme() + ggtitle("Cell Type")
p2 <- DimPlotMetacells(seurat_obj, group.by='Patient') + umap_theme() + ggtitle("Patient")
p1 | p2

saveRDS(seurat_obj, file='hdWGCNA_after-RunUMAPMetacells.rds')

## We previsouly subsetted by celltype, to ensure out data format was consisten with the vignette 
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Astrocytes", # the name of the group of interest in the group.by column
  group.by='Cluster.id' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)

# 7. Select soft-power threshold (identify the power that your correlation matrix should be raised to decrease noise)
## Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # set this to FALSE since we did this above
)
## plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

## assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# 8. Construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=8,
  setDatExpr=FALSE
)
dev.off()

# 9. Module identification
## a 'module" is a group of highly interconnected genes as measured by Topological Overlap Matrix (T.O.M.)
PlotDendrogram(seurat_obj, main='Astrocytes hdWGCNA Dendrogram')

# 10. Compute all Module Eigengenes in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Diagnosis"
)

## Harmonize module eigengenes:
hMEs <- GetMEs(seurat_obj)
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

## Compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Cluster.id', group_name = 'Astrocytes')

## Rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Ast-M"
)

saveRDS(seurat_obj, file='hdWGCNA_Ast_post-eigengenes.rds')

## It is good to generate a separate csv file to keep track of the genes within modules
## get the module assignment table:
modules <- GetModules(seurat_obj)
## show the first 6 columns:
head(modules[,1:6])

write.csv(modules, "Ast_modules.subset-first.csv")


# STEP 2: VISUALIZE THE NETWORK

# 1. Compute gene scoring for the top 25 hub genes by kME for each module
## We used UCell method [https://github.com/carmonalab/UCell]
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# 2. Plot module correlagram
library(igraph)
library(corrplot)

ModuleCorrelogram(seurat_obj)

## add MEvalues to metadata

seurat_obj@meta.data <- cbind(
  seurat_obj@meta.data,
  GetMEs(seurat_obj, harmonized=TRUE)
)

# 2. Plot module eigengene value on "cell levels"
## Plot Ast hME using Seurat VlnPlot function
p1 <- VlnPlot(seurat_obj,
             features = c('Ast-M1'),
             #'Ast-M1','Ast-M2','Ast-M3' 
             #'Ast-M4', 'Ast-M5'
             group.by = 'Diagnosis',
             pt.size = 0)
 p1= p1+geom_boxplot(width=.25, fill='white')

## Change axis labels and remove legend:
p1 <- p1 + xlab('') + ylab('hME') + NoLegend()
## Plot output
p1


# STEP 3. RUN STATISTICAL TESTS 
# See if any particular module is particularly expressed between groups. 

# 1. Build a dataframe of the different module names assigned in network analysis plot.
modules <- c("Ast-M1", "Ast-M2", "Ast-M3", "Ast-M4", "Ast-M5", "Ast-M6", "Ast-M7", "Ast-M8", "Ast-M9", "Ast-M10")
modules <- data.frame(modules)

# A. Kruskal Wallis test
Kruskal_Wallis<- setNames(data.frame(matrix(ncol = 3, nrow = nrow(modules))), c("Module", "Kruskal-Wallis chi-squared", "p-value"))
for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  #Kruskal-Wallis test (one way ANOVA)
  Kruskal_Wallis_data<-kruskal.test(hME ~ diagnosis, data = hME)
  
  Kruskal_Wallis[i,1]<-modules[i,1]
  Kruskal_Wallis[i,2]<-Kruskal_Wallis_data$statistic
  Kruskal_Wallis[i,3]<-Kruskal_Wallis_data$p.value
}
write_csv(Kruskal_Wallis, "Kruskal_Wallis_data.csv")

# B. Wilcoxon rank sum test
## Adapted from: [https://www.datanovia.com/en/lessons/wilcoxon-test-in-r]
library(tidyverse)
library(rstatix)
library(ggpubr)

## Load data. Extract metadata from a particular module you are interested

# PSEN1-E280A vs SporadicAD  
Wilcoxon_rank_sum<- setNames(data.frame(matrix(ncol = 5, nrow = nrow(modules))), c("Module", "Wilcoxon p-value (BH adjusted)", "Delaney's A effect size", "Delaney's A magnitude", "Wilcoxon effect size"))

for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  ## Show a sample of the data by group
  set.seed(223)
  
  ## Signifiance test
  stat.test <- hME %>% 
    rstatix::wilcox_test(hME ~ diagnosis) %>%
    add_significance()
  stat.test
  W_effect_size<-hME %>% wilcox_effsize(hME ~ diagnosis)
  
  Wilcoxon_rank_sum[i,1] <- modules[i,1]
  Wilcoxon_rank_sum[i,2] <- stat.test[1,7]
  Wilcoxon_rank_sum[i,5] <- W_effect_size[1,4]
         
  ## get AUC effect size
  Sporadic<-hME$hME[hME$diagnosis=='Sporadic AD']
  E280A<-hME$hME[hME$diagnosis=='E280A']
         
  Delaneys_A<-VD.A(new_frontal_official, E280A)
  Wilcoxon_rank_sum[i,3] <-Delaneys_A[["estimate"]]
  Wilcoxon_rank_sum[i,4] <-Delaneys_A[["magnitude"]]   
}

write_csv(Wilcoxon_rank_sum, "Wilcoxon_rank_sum_data.csv")

## Graph Wilcoxon rank sum test
for (i in 1:nrow(modules)){
  
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  ## Make boxplot
  bxp <-ggboxplot(
    hME, x = "diagnosis", y = "hME", 
    ylab = "hME", xlab = "Diganosis", add = "jitter"
  )
  ## add to graph
  stat.test <- stat.test %>% add_xy_position(x = "diagnosis")
  bxp <-bxp + 
    stat_pvalue_manual(stat.test, tip.length = 0) +
    labs(title = paste(modules[i,1]))
  assign(paste('bxp_', print(modules[i,1]),sep = ""), bxp)
}

