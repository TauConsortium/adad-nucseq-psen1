#https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/basic_tutorial.Rmd

load("~/new_frontal_official/brainFRONTAL.noMT.RData")

# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)

setwd("~/new_frontal_official/hdWGCNA/Ast-subset")

seurat_obj<-BrainFRONTAL.noMT
seurat_obj<- subset(x = seurat_obj, subset = Cluster.id == "Astrocytes")
DefaultAssay(seurat_obj)<- 'RNA'

seurat_obj<-NormalizeData(seurat_obj)
seurat_obj<-FindVariableFeatures(seurat_obj)
seurat_obj<-ScaleData(seurat_obj)
seurat_obj<-RunPCA(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj<-RunUMAP(seurat_obj, dims=1:20)
seurat_obj<-FindNeighbors(seurat_obj)
seurat_obj<-FindClusters(seurat_obj)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "BrainFP" # the name of the hdWGCNA experiment
)


seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Patient", "Cluster.id"), # specify the columns in seurat_obj@meta.data to group by
  k = 15, # nearest-neighbors parameter, use 15 cells
  ident.group = 'Cluster.id' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)


seurat_obj <- seurat_obj %>%
  NormalizeMetacells() %>%
  ScaleMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunPCAMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunHarmonyMetacells(group.by.vars='Patient') %>%
  RunUMAPMetacells(reduction='harmony', dims=1:15)
p1 <- DimPlotMetacells(seurat_obj, group.by='Cluster.id') + umap_theme() + ggtitle("Cell Type")
p2 <- DimPlotMetacells(seurat_obj, group.by='Patient') + umap_theme() + ggtitle("Patient")
p1 | p2

#setwd("~/new_frontal_official/hdWGCNA-full")

saveRDS(seurat_obj, file='hdWGCNA_after-RunUMAPMetacells.rds')
saveRDS(metacell_obj, file='hdWGCNA_metacell_obj.rds')


#seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/hdWGCNA_after-RunUMAPMetacells.rds")

#Co-expression of specific celltype

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Astrocytes", # the name of the group of interest in the group.by column
  group.by='Cluster.id' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


## Select soft-power threshold
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # set this to FALSE since we did this above
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)



# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=8,
  setDatExpr=FALSE
)

dev.off()

PlotDendrogram(seurat_obj, main='Astrocytes hdWGCNA Dendrogram')


# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Diagnosis"
)


# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Cluster.id', group_name = 'Astrocytes')

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Ast-M"
)

saveRDS(seurat_obj, file='hdWGCNA_Ast_post-eigengenes.rds')

# plot genes ranked by kME for each module
#p <- PlotKMEs(seurat_obj, ncol=5)
#p






#seurat_obj <- readRDS("~/New_Frontal_Pole/hdWGCNA/hdWGCNA_Ast_post-eigengenes.rds")

# get the module assignment table:
modules <- GetModules(seurat_obj)
# show the first 6 columns:
head(modules[,1:6])

write.csv(modules, "Ast_modules.subset-first.csv")

# compute gene scoring for the top 15 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# plot module correlagram
library(igraph)
library(corrplot)

ModuleCorrelogram(seurat_obj)


#add MEvalues to metadata

seurat_obj@meta.data <- cbind(
  seurat_obj@meta.data,
  GetMEs(seurat_obj, harmonized=TRUE)
)

# Plot Ast hME using Seurat VlnPlot function
p1 <- VlnPlot(seurat_obj,
             features = c('Ast-M1'),
             #'Ast-M1','Ast-M2','Ast-M3' 
             #'Ast-M4', 'Ast-M5'
             group.by = 'Diagnosis',
             pt.size = 0)
 p1= p1+geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p1 <- p1 + xlab('') + ylab('hME') + NoLegend()
# plot output
p1


#wilcoxon rank sum test
#https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/#prerequisites
library(tidyverse)
library(rstatix)
library(ggpubr)


hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[["Ast-M1"]])
rownames(hME)=rownames(seurat_obj@meta.data)

# Show a sample of the data by group
set.seed(123)
hME %>% sample_n_by(diagnosis, size = 3)

#get summary data
hME %>%
  group_by(diagnosis) %>%
  get_summary_stats(hME, type = "median_iqr")

#signifiance test
stat.test <- hME %>% 
  rstatix::wilcox_test(hME ~ diagnosis) %>%
  add_significance()
stat.test

hME %>% wilcox_effsize(hME ~ diagnosis)

#make boxplot
library(ggplot2)

bxp <-ggboxplot(
  hME, x = "diagnosis", y = "hME", 
  ylab = "hME", xlab = "Diganosis", add = "jitter"
)
bxp

#add to graph
stat.test <- stat.test %>% add_xy_position(x = "diagnosis")
bxp + 
  stat_pvalue_manual(stat.test, tip.length = 0) +
  labs(title = 'Ast-M1')
get_test_label(stat.test)

saveRDS(seurat_obj, "hdWGCNA-seurat_obj-basic-end.rds")



seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/Ast-subset/hdWGCNA-seurat_obj-basic-end.rds")
