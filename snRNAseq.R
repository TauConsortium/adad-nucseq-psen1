library(dplyr)
library(Seurat)
library (Matrix)
library(ggplot2)
library(sctransform)
library (EnhancedVolcano)
library (DoubletFinder)
library (pheatmap)

####Pre-processing ----
####Load Dataset
C114_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C114_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C114_C01 <- CreateSeuratObject(counts = C114_C01.data, project = "C114_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C114_C01 <- PercentageFeatureSet(C114_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C114_C01)
VlnPlot(C114_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C114_C01 <- subset(C114_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C114_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C114_C01)
## Pre-process Seurat object (sctransform)
C114_C01 <- SCTransform(C114_C01)
C114_C01 <- RunPCA(C114_C01)
C114_C01 <- RunUMAP(C114_C01, dims = 1:20)
####Doublet Finder----
##114_C01
## pK Identification (no ground-truth)
sweep.res.list_C114_C01 <- paramSweep_v3(C114_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C114_C01 <- summarizeSweep(sweep.res.list_C114_C01, GT = FALSE)
bcmvn_C114_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C114_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C114_C01 <- doubletFinder_v3(C114_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C114_C01 <- doubletFinder_v3(C114_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C114_C01@meta.data[,"DF_hi.lo"] <- C114_C01@meta.data$DF.classifications_0.25_0.09_82
C114_C01@meta.data$DF_hi.lo[which(C114_C01@meta.data$DF_hi.lo == "Doublet" & C114_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C114_C01@meta.data$DF_hi.lo[which(C114_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets----
#114_C01
dim (C114_C01)
Idents(C114_C01) <- "DF_hi.lo"
C114_C01 <- subset(C114_C01, idents = "Singlet")
dim (C114_C01)
head (C114_C01@meta.data)
tail (C114_C01@meta.data)
####Remove unwanted genes (MT-)----
#C114_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C114_C01)
total.genes <- list(rownames(C114_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C114_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C114_C01 <- subset (C114_C01, features = total.genes_subset)
dim(C114_C01)

####Pre-processing ----
####Load Dataset
C135_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C135_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C135_C01 <- CreateSeuratObject(counts = C135_C01.data, project = "C135_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C135_C01 <- PercentageFeatureSet(C135_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C135_C01)
VlnPlot(C135_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C135_C01 <- subset(C135_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C135_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C135_C01)
## Pre-process Seurat object (sctransform)
C135_C01 <- SCTransform(C135_C01)
C135_C01 <- RunPCA(C135_C01)
C135_C01 <- RunUMAP(C135_C01, dims = 1:20)
####Doublet Finder
##C135_C01
## pK Identification (no ground-truth)
sweep.res.list_C135_C01 <- paramSweep_v3(C135_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C135_C01 <- summarizeSweep(sweep.res.list_C135_C01, GT = FALSE)
bcmvn_C135_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C135_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C135_C01 <- doubletFinder_v3(C135_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C135_C01 <- doubletFinder_v3(C135_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C135_C01@meta.data[,"DF_hi.lo"] <- C135_C01@meta.data$DF.classifications_0.25_0.09_82
C135_C01@meta.data$DF_hi.lo[which(C135_C01@meta.data$DF_hi.lo == "Doublet" & C135_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C135_C01@meta.data$DF_hi.lo[which(C135_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#135_C01
dim (C135_C01)
Idents(C135_C01) <- "DF_hi.lo"
C135_C01 <- subset(C135_C01, idents = "Singlet")
dim (C135_C01)
head (C135_C01@meta.data)
tail (C135_C01@meta.data)
####Remove unwanted genes (MT-)
#C135_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C135_C01)
total.genes <- list(rownames(C135_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C135_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C135_C01 <- subset (C135_C01, features = total.genes_subset)
dim(C135_C01)

####Pre-processing ----
####Load Dataset
C140_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C140_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C140_C01 <- CreateSeuratObject(counts = C140_C01.data, project = "C140_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C140_C01 <- PercentageFeatureSet(C140_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C140_C01)
VlnPlot(C140_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C140_C01 <- subset(C140_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C140_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C140_C01)
## Pre-process Seurat object (sctransform)
C140_C01 <- SCTransform(C140_C01)
C140_C01 <- RunPCA(C140_C01)
C140_C01 <- RunUMAP(C140_C01, dims = 1:20)
####Doublet Finder
##C140_C01
## pK Identification (no ground-truth)
sweep.res.list_C140_C01 <- paramSweep_v3(C140_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C140_C01 <- summarizeSweep(sweep.res.list_C140_C01, GT = FALSE)
bcmvn_C140_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C140_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C140_C01 <- doubletFinder_v3(C140_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C140_C01 <- doubletFinder_v3(C140_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C140_C01@meta.data[,"DF_hi.lo"] <- C140_C01@meta.data$DF.classifications_0.25_0.09_82
C140_C01@meta.data$DF_hi.lo[which(C140_C01@meta.data$DF_hi.lo == "Doublet" & C140_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C140_C01@meta.data$DF_hi.lo[which(C140_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#140_C01
dim (C140_C01)
Idents(C140_C01) <- "DF_hi.lo"
C140_C01 <- subset(C140_C01, idents = "Singlet")
dim (C140_C01)
head (C140_C01@meta.data)
tail (C140_C01@meta.data)
####Remove unwanted genes (MT-)
#C140_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C140_C01)
total.genes <- list(rownames(C140_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C140_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C140_C01 <- subset (C140_C01, features = total.genes_subset)
dim(C140_C01)

####Pre-processing ----
####Load Dataset
C142_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C142_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C142_C01 <- CreateSeuratObject(counts = C142_C01.data, project = "C142_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C142_C01 <- PercentageFeatureSet(C142_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C142_C01)
VlnPlot(C142_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C142_C01 <- subset(C142_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C142_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C142_C01)
## Pre-process Seurat object (sctransform)
C142_C01 <- SCTransform(C142_C01)
C142_C01 <- RunPCA(C142_C01)
C142_C01 <- RunUMAP(C142_C01, dims = 1:20)
####Doublet Finder
##C142_C01
## pK Identification (no ground-truth)
sweep.res.list_C142_C01 <- paramSweep_v3(C142_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C142_C01 <- summarizeSweep(sweep.res.list_C142_C01, GT = FALSE)
bcmvn_C142_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C142_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C142_C01 <- doubletFinder_v3(C142_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C142_C01 <- doubletFinder_v3(C142_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C142_C01@meta.data[,"DF_hi.lo"] <- C142_C01@meta.data$DF.classifications_0.25_0.09_82
C142_C01@meta.data$DF_hi.lo[which(C142_C01@meta.data$DF_hi.lo == "Doublet" & C142_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C142_C01@meta.data$DF_hi.lo[which(C142_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#142_C01
dim (C142_C01)
Idents(C142_C01) <- "DF_hi.lo"
C142_C01 <- subset(C142_C01, idents = "Singlet")
dim (C142_C01)
head (C142_C01@meta.data)
tail (C142_C01@meta.data)
####Remove unwanted genes (MT-)
#C142_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C142_C01)
total.genes <- list(rownames(C142_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C142_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C142_C01 <- subset (C142_C01, features = total.genes_subset)
dim(C142_C01)

####Pre-processing ----
####Load Dataset
C361_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C361_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C361_C01 <- CreateSeuratObject(counts = C361_C01.data, project = "C361_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C361_C01 <- PercentageFeatureSet(C361_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C361_C01)
VlnPlot(C361_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C361_C01 <- subset(C361_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C361_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C361_C01)
## Pre-process Seurat object (sctransform)
C361_C01 <- SCTransform(C361_C01)
C361_C01 <- RunPCA(C361_C01)
C361_C01 <- RunUMAP(C361_C01, dims = 1:20)
####Doublet Finder
##C361_C01
## pK Identification (no ground-truth)
sweep.res.list_C361_C01 <- paramSweep_v3(C361_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C361_C01 <- summarizeSweep(sweep.res.list_C361_C01, GT = FALSE)
bcmvn_C361_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C361_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C361_C01 <- doubletFinder_v3(C361_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C361_C01 <- doubletFinder_v3(C361_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C361_C01@meta.data[,"DF_hi.lo"] <- C361_C01@meta.data$DF.classifications_0.25_0.09_82
C361_C01@meta.data$DF_hi.lo[which(C361_C01@meta.data$DF_hi.lo == "Doublet" & C361_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C361_C01@meta.data$DF_hi.lo[which(C361_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#361_C01
dim (C361_C01)
Idents(C361_C01) <- "DF_hi.lo"
C361_C01 <- subset(C361_C01, idents = "Singlet")
dim (C361_C01)
head (C361_C01@meta.data)
tail (C361_C01@meta.data)
####Remove unwanted genes (MT-)
#C361_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C361_C01)
total.genes <- list(rownames(C361_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C361_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C361_C01 <- subset (C361_C01, features = total.genes_subset)
dim(C361_C01)

####Pre-processing ----
####Load Dataset
C363_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C363_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C363_C01 <- CreateSeuratObject(counts = C363_C01.data, project = "C363_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C363_C01 <- PercentageFeatureSet(C363_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C363_C01)
VlnPlot(C363_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C363_C01 <- subset(C363_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C363_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C363_C01)
## Pre-process Seurat object (sctransform)
C363_C01 <- SCTransform(C363_C01)
C363_C01 <- RunPCA(C363_C01)
C363_C01 <- RunUMAP(C363_C01, dims = 1:20)
####Doublet Finder
##C363_C01
## pK Identification (no ground-truth)
sweep.res.list_C363_C01 <- paramSweep_v3(C363_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C363_C01 <- summarizeSweep(sweep.res.list_C363_C01, GT = FALSE)
bcmvn_C363_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C363_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C363_C01 <- doubletFinder_v3(C363_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C363_C01 <- doubletFinder_v3(C363_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C363_C01@meta.data[,"DF_hi.lo"] <- C363_C01@meta.data$DF.classifications_0.25_0.09_82
C363_C01@meta.data$DF_hi.lo[which(C363_C01@meta.data$DF_hi.lo == "Doublet" & C363_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C363_C01@meta.data$DF_hi.lo[which(C363_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#363_C01
dim (C363_C01)
Idents(C363_C01) <- "DF_hi.lo"
C363_C01 <- subset(C363_C01, idents = "Singlet")
dim (C363_C01)
head (C363_C01@meta.data)
tail (C363_C01@meta.data)
####Remove unwanted genes (MT-)
#C363_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C363_C01)
total.genes <- list(rownames(C363_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C363_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C363_C01 <- subset (C363_C01, features = total.genes_subset)
dim(C363_C01)

####Pre-processing ----
####Load Dataset
C379_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C379_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C379_C01 <- CreateSeuratObject(counts = C379_C01.data, project = "C379_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C379_C01 <- PercentageFeatureSet(C379_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C379_C01)
VlnPlot(C379_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C379_C01 <- subset(C379_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C379_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C379_C01)
## Pre-process Seurat object (sctransform)
C379_C01 <- SCTransform(C379_C01)
C379_C01 <- RunPCA(C379_C01)
C379_C01 <- RunUMAP(C379_C01, dims = 1:20)
####Doublet Finder
##C379_C01
## pK Identification (no ground-truth)
sweep.res.list_C379_C01 <- paramSweep_v3(C379_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C379_C01 <- summarizeSweep(sweep.res.list_C379_C01, GT = FALSE)
bcmvn_C379_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C379_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C379_C01 <- doubletFinder_v3(C379_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C379_C01 <- doubletFinder_v3(C379_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C379_C01@meta.data[,"DF_hi.lo"] <- C379_C01@meta.data$DF.classifications_0.25_0.09_82
C379_C01@meta.data$DF_hi.lo[which(C379_C01@meta.data$DF_hi.lo == "Doublet" & C379_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C379_C01@meta.data$DF_hi.lo[which(C379_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#379_C01
dim (C379_C01)
Idents(C379_C01) <- "DF_hi.lo"
C379_C01 <- subset(C379_C01, idents = "Singlet")
dim (C379_C01)
head (C379_C01@meta.data)
tail (C379_C01@meta.data)
####Remove unwanted genes (MT-)
#C379_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C379_C01)
total.genes <- list(rownames(C379_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C379_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C379_C01 <- subset (C379_C01, features = total.genes_subset)
dim(C379_C01)

####Pre-processing ----
####Load Dataset
C377_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C377_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C377_C01 <- CreateSeuratObject(counts = C377_C01.data, project = "C377_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C377_C01 <- PercentageFeatureSet(C377_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C377_C01)
VlnPlot(C377_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C377_C01 <- subset(C377_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C377_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C377_C01)
## Pre-process Seurat object (sctransform)
C377_C01 <- SCTransform(C377_C01)
C377_C01 <- RunPCA(C377_C01)
C377_C01 <- RunUMAP(C377_C01, dims = 1:20)
####Doublet Finder
##C377_C01
## pK Identification (no ground-truth)
sweep.res.list_C377_C01 <- paramSweep_v3(C377_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C377_C01 <- summarizeSweep(sweep.res.list_C377_C01, GT = FALSE)
bcmvn_C377_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C377_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C377_C01 <- doubletFinder_v3(C377_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C377_C01 <- doubletFinder_v3(C377_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C377_C01@meta.data[,"DF_hi.lo"] <- C377_C01@meta.data$DF.classifications_0.25_0.09_82
C377_C01@meta.data$DF_hi.lo[which(C377_C01@meta.data$DF_hi.lo == "Doublet" & C377_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C377_C01@meta.data$DF_hi.lo[which(C377_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#377_C01
dim (C377_C01)
Idents(C377_C01) <- "DF_hi.lo"
C377_C01 <- subset(C377_C01, idents = "Singlet")
dim (C377_C01)
head (C377_C01@meta.data)
tail (C377_C01@meta.data)
####Remove unwanted genes (MT-)
#C377_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C377_C01)
total.genes <- list(rownames(C377_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C377_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C377_C01 <- subset (C377_C01, features = total.genes_subset)
dim(C377_C01)

####Pre-processing ----
####Load Dataset
C381_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C381_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C381_C01 <- CreateSeuratObject(counts = C381_C01.data, project = "C381_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C381_C01 <- PercentageFeatureSet(C381_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C381_C01)
VlnPlot(C381_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C381_C01 <- subset(C381_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C381_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C381_C01)
## Pre-process Seurat object (sctransform)
C381_C01 <- SCTransform(C381_C01)
C381_C01 <- RunPCA(C381_C01)
C381_C01 <- RunUMAP(C381_C01, dims = 1:20)
####Doublet Finder
##C381_C01
## pK Identification (no ground-truth)
sweep.res.list_C381_C01 <- paramSweep_v3(C381_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C381_C01 <- summarizeSweep(sweep.res.list_C381_C01, GT = FALSE)
bcmvn_C381_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C381_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C381_C01 <- doubletFinder_v3(C381_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C381_C01 <- doubletFinder_v3(C381_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C381_C01@meta.data[,"DF_hi.lo"] <- C381_C01@meta.data$DF.classifications_0.25_0.09_82
C381_C01@meta.data$DF_hi.lo[which(C381_C01@meta.data$DF_hi.lo == "Doublet" & C381_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C381_C01@meta.data$DF_hi.lo[which(C381_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#381_C01
dim (C381_C01)
Idents(C381_C01) <- "DF_hi.lo"
C381_C01 <- subset(C381_C01, idents = "Singlet")
dim (C381_C01)
head (C381_C01@meta.data)
tail (C381_C01@meta.data)
####Remove unwanted genes (MT-)
#C381_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C381_C01)
total.genes <- list(rownames(C381_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C381_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C381_C01 <- subset (C381_C01, features = total.genes_subset)
dim(C381_C01)

####Pre-processing ----
####Load Dataset
C323_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C323_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C323_C01 <- CreateSeuratObject(counts = C323_C01.data, project = "C323_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C323_C01 <- PercentageFeatureSet(C323_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C323_C01)
VlnPlot(C323_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C323_C01 <- subset(C323_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C323_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C323_C01)
## Pre-process Seurat object (sctransform)
C323_C01 <- SCTransform(C323_C01)
C323_C01 <- RunPCA(C323_C01)
C323_C01 <- RunUMAP(C323_C01, dims = 1:20)
####Doublet Finder
##C323_C01
## pK Identification (no ground-truth)
sweep.res.list_C323_C01 <- paramSweep_v3(C323_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C323_C01 <- summarizeSweep(sweep.res.list_C323_C01, GT = FALSE)
bcmvn_C323_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C323_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C323_C01 <- doubletFinder_v3(C323_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C323_C01 <- doubletFinder_v3(C323_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C323_C01@meta.data[,"DF_hi.lo"] <- C323_C01@meta.data$DF.classifications_0.25_0.09_82
C323_C01@meta.data$DF_hi.lo[which(C323_C01@meta.data$DF_hi.lo == "Doublet" & C323_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C323_C01@meta.data$DF_hi.lo[which(C323_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#323_C01
dim (C323_C01)
Idents(C323_C01) <- "DF_hi.lo"
C323_C01 <- subset(C323_C01, idents = "Singlet")
dim (C323_C01)
head (C323_C01@meta.data)
tail (C323_C01@meta.data)
####Remove unwanted genes (MT-)
#C323_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C323_C01)
total.genes <- list(rownames(C323_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C323_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C323_C01 <- subset (C323_C01, features = total.genes_subset)
dim(C323_C01)

####Pre-processing ----
####Load Dataset
C255_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C255_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C255_C01 <- CreateSeuratObject(counts = C255_C01.data, project = "C255_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C255_C01 <- PercentageFeatureSet(C255_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C255_C01)
VlnPlot(C255_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C255_C01 <- subset(C255_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C255_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C255_C01)
## Pre-process Seurat object (sctransform)
C255_C01 <- SCTransform(C255_C01)
C255_C01 <- RunPCA(C255_C01)
C255_C01 <- RunUMAP(C255_C01, dims = 1:20)
####Doublet Finder
##C255_C01
## pK Identification (no ground-truth)
sweep.res.list_C255_C01 <- paramSweep_v3(C255_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C255_C01 <- summarizeSweep(sweep.res.list_C255_C01, GT = FALSE)
bcmvn_C255_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C255_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C255_C01 <- doubletFinder_v3(C255_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C255_C01 <- doubletFinder_v3(C255_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C255_C01@meta.data[,"DF_hi.lo"] <- C255_C01@meta.data$DF.classifications_0.25_0.09_82
C255_C01@meta.data$DF_hi.lo[which(C255_C01@meta.data$DF_hi.lo == "Doublet" & C255_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C255_C01@meta.data$DF_hi.lo[which(C255_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#255_C01
dim (C255_C01)
Idents(C255_C01) <- "DF_hi.lo"
C255_C01 <- subset(C255_C01, idents = "Singlet")
dim (C255_C01)
head (C255_C01@meta.data)
tail (C255_C01@meta.data)
####Remove unwanted genes (MT-)
#C255_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C255_C01)
total.genes <- list(rownames(C255_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C255_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C255_C01 <- subset (C255_C01, features = total.genes_subset)
dim(C255_C01)

####Pre-processing ----
####Load Dataset
C210_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C210_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C210_C01 <- CreateSeuratObject(counts = C210_C01.data, project = "C210_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C210_C01 <- PercentageFeatureSet(C210_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C210_C01)
VlnPlot(C210_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C210_C01 <- subset(C210_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C210_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C210_C01)
## Pre-process Seurat object (sctransform)
C210_C01 <- SCTransform(C210_C01)
C210_C01 <- RunPCA(C210_C01)
C210_C01 <- RunUMAP(C210_C01, dims = 1:20)
####Doublet Finder
##C210_C01
## pK Identification (no ground-truth)
sweep.res.list_C210_C01 <- paramSweep_v3(C210_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C210_C01 <- summarizeSweep(sweep.res.list_C210_C01, GT = FALSE)
bcmvn_C210_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C210_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C210_C01 <- doubletFinder_v3(C210_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C210_C01 <- doubletFinder_v3(C210_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C210_C01@meta.data[,"DF_hi.lo"] <- C210_C01@meta.data$DF.classifications_0.25_0.09_82
C210_C01@meta.data$DF_hi.lo[which(C210_C01@meta.data$DF_hi.lo == "Doublet" & C210_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C210_C01@meta.data$DF_hi.lo[which(C210_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#210_C01
dim (C210_C01)
Idents(C210_C01) <- "DF_hi.lo"
C210_C01 <- subset(C210_C01, idents = "Singlet")
dim (C210_C01)
head (C210_C01@meta.data)
tail (C210_C01@meta.data)
####Remove unwanted genes (MT-)
#C210_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C210_C01)
total.genes <- list(rownames(C210_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C210_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C210_C01 <- subset (C210_C01, features = total.genes_subset)
dim(C210_C01)

####Pre-processing ----
####Load Dataset
C269_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C269_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C269_C01 <- CreateSeuratObject(counts = C269_C01.data, project = "C269_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C269_C01 <- PercentageFeatureSet(C269_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C269_C01)
VlnPlot(C269_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C269_C01 <- subset(C269_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C269_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C269_C01)
## Pre-process Seurat object (sctransform)
C269_C01 <- SCTransform(C269_C01)
C269_C01 <- RunPCA(C269_C01)
C269_C01 <- RunUMAP(C269_C01, dims = 1:20)
####Doublet Finder
##C269_C01
## pK Identification (no ground-truth)
sweep.res.list_C269_C01 <- paramSweep_v3(C269_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C269_C01 <- summarizeSweep(sweep.res.list_C269_C01, GT = FALSE)
bcmvn_C269_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C269_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C269_C01 <- doubletFinder_v3(C269_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C269_C01 <- doubletFinder_v3(C269_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C269_C01@meta.data[,"DF_hi.lo"] <- C269_C01@meta.data$DF.classifications_0.25_0.09_82
C269_C01@meta.data$DF_hi.lo[which(C269_C01@meta.data$DF_hi.lo == "Doublet" & C269_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C269_C01@meta.data$DF_hi.lo[which(C269_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#269_C01
dim (C269_C01)
Idents(C269_C01) <- "DF_hi.lo"
C269_C01 <- subset(C269_C01, idents = "Singlet")
dim (C269_C01)
head (C269_C01@meta.data)
tail (C269_C01@meta.data)
####Remove unwanted genes (MT-)
#C269_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C269_C01)
total.genes <- list(rownames(C269_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C269_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C269_C01 <- subset (C269_C01, features = total.genes_subset)
dim(C269_C01)

####Pre-processing ----
####Load Dataset
C99_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C99_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C99_C01 <- CreateSeuratObject(counts = C99_C01.data, project = "C99_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C99_C01 <- PercentageFeatureSet(C99_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C99_C01)
VlnPlot(C99_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C99_C01 <- subset(C99_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C99_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C99_C01)
## Pre-process Seurat object (sctransform)
C99_C01 <- SCTransform(C99_C01)
C99_C01 <- RunPCA(C99_C01)
C99_C01 <- RunUMAP(C99_C01, dims = 1:20)
####Doublet Finder
##C99_C01
## pK Identification (no ground-truth)
sweep.res.list_C99_C01 <- paramSweep_v3(C99_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C99_C01 <- summarizeSweep(sweep.res.list_C99_C01, GT = FALSE)
bcmvn_C99_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C99_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C99_C01 <- doubletFinder_v3(C99_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C99_C01 <- doubletFinder_v3(C99_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C99_C01@meta.data[,"DF_hi.lo"] <- C99_C01@meta.data$DF.classifications_0.25_0.09_82
C99_C01@meta.data$DF_hi.lo[which(C99_C01@meta.data$DF_hi.lo == "Doublet" & C99_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C99_C01@meta.data$DF_hi.lo[which(C99_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#99_C01
dim (C99_C01)
Idents(C99_C01) <- "DF_hi.lo"
C99_C01 <- subset(C99_C01, idents = "Singlet")
dim (C99_C01)
head (C99_C01@meta.data)
tail (C99_C01@meta.data)
####Remove unwanted genes (MT-)
#C99_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C99_C01)
total.genes <- list(rownames(C99_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C99_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C99_C01 <- subset (C99_C01, features = total.genes_subset)
dim(C99_C01)

####Pre-processing ----
####Load Dataset
C181_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C181_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C181_C01 <- CreateSeuratObject(counts = C181_C01.data, project = "C181_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C181_C01 <- PercentageFeatureSet(C181_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C181_C01)
VlnPlot(C181_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C181_C01 <- subset(C181_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C181_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C181_C01)
## Pre-process Seurat object (sctransform)
C181_C01 <- SCTransform(C181_C01)
C181_C01 <- RunPCA(C181_C01)
C181_C01 <- RunUMAP(C181_C01, dims = 1:20)
####Doublet Finder
##C181_C01
## pK Identification (no ground-truth)
sweep.res.list_C181_C01 <- paramSweep_v3(C181_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C181_C01 <- summarizeSweep(sweep.res.list_C181_C01, GT = FALSE)
bcmvn_C181_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C181_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C181_C01 <- doubletFinder_v3(C181_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C181_C01 <- doubletFinder_v3(C181_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C181_C01@meta.data[,"DF_hi.lo"] <- C181_C01@meta.data$DF.classifications_0.25_0.09_82
C181_C01@meta.data$DF_hi.lo[which(C181_C01@meta.data$DF_hi.lo == "Doublet" & C181_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C181_C01@meta.data$DF_hi.lo[which(C181_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#181_C01
dim (C181_C01)
Idents(C181_C01) <- "DF_hi.lo"
C181_C01 <- subset(C181_C01, idents = "Singlet")
dim (C181_C01)
head (C181_C01@meta.data)
tail (C181_C01@meta.data)
####Remove unwanted genes (MT-)
#C181_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C181_C01)
total.genes <- list(rownames(C181_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C181_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C181_C01 <- subset (C181_C01, features = total.genes_subset)
dim(C181_C01)

####Pre-processing ----
####Load Dataset
C277_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C277_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C277_C01 <- CreateSeuratObject(counts = C277_C01.data, project = "C277_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C277_C01 <- PercentageFeatureSet(C277_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C277_C01)
VlnPlot(C277_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C277_C01 <- subset(C277_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C277_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C277_C01)
## Pre-process Seurat object (sctransform)
C277_C01 <- SCTransform(C277_C01)
C277_C01 <- RunPCA(C277_C01)
C277_C01 <- RunUMAP(C277_C01, dims = 1:20)
####Doublet Finder
##C277_C01
## pK Identification (no ground-truth)
sweep.res.list_C277_C01 <- paramSweep_v3(C277_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C277_C01 <- summarizeSweep(sweep.res.list_C277_C01, GT = FALSE)
bcmvn_C277_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C277_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C277_C01 <- doubletFinder_v3(C277_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C277_C01 <- doubletFinder_v3(C277_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C277_C01@meta.data[,"DF_hi.lo"] <- C277_C01@meta.data$DF.classifications_0.25_0.09_82
C277_C01@meta.data$DF_hi.lo[which(C277_C01@meta.data$DF_hi.lo == "Doublet" & C277_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C277_C01@meta.data$DF_hi.lo[which(C277_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#277_C01
dim (C277_C01)
Idents(C277_C01) <- "DF_hi.lo"
C277_C01 <- subset(C277_C01, idents = "Singlet")
dim (C277_C01)
head (C277_C01@meta.data)
tail (C277_C01@meta.data)
####Remove unwanted genes (MT-)
#C277_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C277_C01)
total.genes <- list(rownames(C277_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C277_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C277_C01 <- subset (C277_C01, features = total.genes_subset)
dim(C277_C01)

####Pre-processing ----
####Load Dataset
C226_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C226_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C226_C01 <- CreateSeuratObject(counts = C226_C01.data, project = "C226_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C226_C01 <- PercentageFeatureSet(C226_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C226_C01)
VlnPlot(C226_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C226_C01 <- subset(C226_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C226_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C226_C01)
## Pre-process Seurat object (sctransform)
C226_C01 <- SCTransform(C226_C01)
C226_C01 <- RunPCA(C226_C01)
C226_C01 <- RunUMAP(C226_C01, dims = 1:20)
####Doublet Finder
##C226_C01
## pK Identification (no ground-truth)
sweep.res.list_C226_C01 <- paramSweep_v3(C226_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C226_C01 <- summarizeSweep(sweep.res.list_C226_C01, GT = FALSE)
bcmvn_C226_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C226_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C226_C01 <- doubletFinder_v3(C226_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C226_C01 <- doubletFinder_v3(C226_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C226_C01@meta.data[,"DF_hi.lo"] <- C226_C01@meta.data$DF.classifications_0.25_0.09_82
C226_C01@meta.data$DF_hi.lo[which(C226_C01@meta.data$DF_hi.lo == "Doublet" & C226_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C226_C01@meta.data$DF_hi.lo[which(C226_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#226_C01
dim (C226_C01)
Idents(C226_C01) <- "DF_hi.lo"
C226_C01 <- subset(C226_C01, idents = "Singlet")
dim (C226_C01)
head (C226_C01@meta.data)
tail (C226_C01@meta.data)
####Remove unwanted genes (MT-)
#C226_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C226_C01)
total.genes <- list(rownames(C226_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C226_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C226_C01 <- subset (C226_C01, features = total.genes_subset)
dim(C226_C01)

####Pre-processing ----
####Load Dataset
C240_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C240_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C240_C01 <- CreateSeuratObject(counts = C240_C01.data, project = "C240_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C240_C01 <- PercentageFeatureSet(C240_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C240_C01)
VlnPlot(C240_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C240_C01 <- subset(C240_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C240_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C240_C01)
## Pre-process Seurat object (sctransform)
C240_C01 <- SCTransform(C240_C01)
C240_C01 <- RunPCA(C240_C01)
C240_C01 <- RunUMAP(C240_C01, dims = 1:20)
####Doublet Finder
##C240_C01
## pK Identification (no ground-truth)
sweep.res.list_C240_C01 <- paramSweep_v3(C240_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C240_C01 <- summarizeSweep(sweep.res.list_C240_C01, GT = FALSE)
bcmvn_C240_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C240_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C240_C01 <- doubletFinder_v3(C240_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C240_C01 <- doubletFinder_v3(C240_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C240_C01@meta.data[,"DF_hi.lo"] <- C240_C01@meta.data$DF.classifications_0.25_0.09_82
C240_C01@meta.data$DF_hi.lo[which(C240_C01@meta.data$DF_hi.lo == "Doublet" & C240_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C240_C01@meta.data$DF_hi.lo[which(C240_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#240_C01
dim (C240_C01)
Idents(C240_C01) <- "DF_hi.lo"
C240_C01 <- subset(C240_C01, idents = "Singlet")
dim (C240_C01)
head (C240_C01@meta.data)
tail (C240_C01@meta.data)
####Remove unwanted genes (MT-)
#C240_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C240_C01)
total.genes <- list(rownames(C240_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C240_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C240_C01 <- subset (C240_C01, features = total.genes_subset)
dim(C240_C01)

####Pre-processing ----
####Load Dataset
C262_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C262_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C262_C01 <- CreateSeuratObject(counts = C262_C01.data, project = "C262_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C262_C01 <- PercentageFeatureSet(C262_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C262_C01)
VlnPlot(C262_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C262_C01 <- subset(C262_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C262_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C262_C01)
## Pre-process Seurat object (sctransform)
C262_C01 <- SCTransform(C262_C01)
C262_C01 <- RunPCA(C262_C01)
C262_C01 <- RunUMAP(C262_C01, dims = 1:20)
####Doublet Finder
##C262_C01
## pK Identification (no ground-truth)
sweep.res.list_C262_C01 <- paramSweep_v3(C262_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C262_C01 <- summarizeSweep(sweep.res.list_C262_C01, GT = FALSE)
bcmvn_C262_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C262_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C262_C01 <- doubletFinder_v3(C262_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C262_C01 <- doubletFinder_v3(C262_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C262_C01@meta.data[,"DF_hi.lo"] <- C262_C01@meta.data$DF.classifications_0.25_0.09_82
C262_C01@meta.data$DF_hi.lo[which(C262_C01@meta.data$DF_hi.lo == "Doublet" & C262_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C262_C01@meta.data$DF_hi.lo[which(C262_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#262_C01
dim (C262_C01)
Idents(C262_C01) <- "DF_hi.lo"
C262_C01 <- subset(C262_C01, idents = "Singlet")
dim (C262_C01)
head (C262_C01@meta.data)
tail (C262_C01@meta.data)
####Remove unwanted genes (MT-)
#C262_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C262_C01)
total.genes <- list(rownames(C262_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C262_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C262_C01 <- subset (C262_C01, features = total.genes_subset)
dim(C262_C01)

####Pre-processing ----
####Load Dataset
C291_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C291_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C291_C01 <- CreateSeuratObject(counts = C291_C01.data, project = "C291_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C291_C01 <- PercentageFeatureSet(C291_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C291_C01)
VlnPlot(C291_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C291_C01 <- subset(C291_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C291_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C291_C01)
## Pre-process Seurat object (sctransform)
C291_C01 <- SCTransform(C291_C01)
C291_C01 <- RunPCA(C291_C01)
C291_C01 <- RunUMAP(C291_C01, dims = 1:20)
####Doublet Finder
##C291_C01
## pK Identification (no ground-truth)
sweep.res.list_C291_C01 <- paramSweep_v3(C291_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C291_C01 <- summarizeSweep(sweep.res.list_C291_C01, GT = FALSE)
bcmvn_C291_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C291_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C291_C01 <- doubletFinder_v3(C291_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C291_C01 <- doubletFinder_v3(C291_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C291_C01@meta.data[,"DF_hi.lo"] <- C291_C01@meta.data$DF.classifications_0.25_0.09_82
C291_C01@meta.data$DF_hi.lo[which(C291_C01@meta.data$DF_hi.lo == "Doublet" & C291_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C291_C01@meta.data$DF_hi.lo[which(C291_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#291_C01
dim (C291_C01)
Idents(C291_C01) <- "DF_hi.lo"
C291_C01 <- subset(C291_C01, idents = "Singlet")
dim (C291_C01)
head (C291_C01@meta.data)
tail (C291_C01@meta.data)
####Remove unwanted genes (MT-)
#C291_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C291_C01)
total.genes <- list(rownames(C291_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C291_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C291_C01 <- subset (C291_C01, features = total.genes_subset)
dim(C291_C01)

####Pre-processing ----
####Load Dataset
C327_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C327_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C327_C01 <- CreateSeuratObject(counts = C327_C01.data, project = "C327_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C327_C01 <- PercentageFeatureSet(C327_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C327_C01)
VlnPlot(C327_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C327_C01 <- subset(C327_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C327_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C327_C01)
## Pre-process Seurat object (sctransform)
C327_C01 <- SCTransform(C327_C01)
C327_C01 <- RunPCA(C327_C01)
C327_C01 <- RunUMAP(C327_C01, dims = 1:20)
####Doublet Finder
##C327_C01
## pK Identification (no ground-truth)
sweep.res.list_C327_C01 <- paramSweep_v3(C327_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C327_C01 <- summarizeSweep(sweep.res.list_C327_C01, GT = FALSE)
bcmvn_C327_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C327_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C327_C01 <- doubletFinder_v3(C327_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C327_C01 <- doubletFinder_v3(C327_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C327_C01@meta.data[,"DF_hi.lo"] <- C327_C01@meta.data$DF.classifications_0.25_0.09_82
C327_C01@meta.data$DF_hi.lo[which(C327_C01@meta.data$DF_hi.lo == "Doublet" & C327_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C327_C01@meta.data$DF_hi.lo[which(C327_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#327_C01
dim (C327_C01)
Idents(C327_C01) <- "DF_hi.lo"
C327_C01 <- subset(C327_C01, idents = "Singlet")
dim (C327_C01)
head (C327_C01@meta.data)
tail (C327_C01@meta.data)
####Remove unwanted genes (MT-)
#C327_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C327_C01)
total.genes <- list(rownames(C327_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C327_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C327_C01 <- subset (C327_C01, features = total.genes_subset)
dim(C327_C01)

####Pre-processing ----
####Load Dataset
C342_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C342_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C342_C01 <- CreateSeuratObject(counts = C342_C01.data, project = "C342_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C342_C01 <- PercentageFeatureSet(C342_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C342_C01)
VlnPlot(C342_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C342_C01 <- subset(C342_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C342_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C342_C01)
## Pre-process Seurat object (sctransform)
C342_C01 <- SCTransform(C342_C01)
C342_C01 <- RunPCA(C342_C01)
C342_C01 <- RunUMAP(C342_C01, dims = 1:20)
####Doublet Finder
##C342_C01
## pK Identification (no ground-truth)
sweep.res.list_C342_C01 <- paramSweep_v3(C342_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C342_C01 <- summarizeSweep(sweep.res.list_C342_C01, GT = FALSE)
bcmvn_C342_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C342_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C342_C01 <- doubletFinder_v3(C342_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C342_C01 <- doubletFinder_v3(C342_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C342_C01@meta.data[,"DF_hi.lo"] <- C342_C01@meta.data$DF.classifications_0.25_0.09_82
C342_C01@meta.data$DF_hi.lo[which(C342_C01@meta.data$DF_hi.lo == "Doublet" & C342_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C342_C01@meta.data$DF_hi.lo[which(C342_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#342_C01
dim (C342_C01)
Idents(C342_C01) <- "DF_hi.lo"
C342_C01 <- subset(C342_C01, idents = "Singlet")
dim (C342_C01)
head (C342_C01@meta.data)
tail (C342_C01@meta.data)
####Remove unwanted genes (MT-)
#C342_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C342_C01)
total.genes <- list(rownames(C342_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C342_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C342_C01 <- subset (C342_C01, features = total.genes_subset)
dim(C342_C01)

####Pre-processing ----
####Load Dataset
C350_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C350_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C350_C01 <- CreateSeuratObject(counts = C350_C01.data, project = "C350_C01", min.cells = 0, min.features = 0)
#### Store mt.percent
C350_C01 <- PercentageFeatureSet(C350_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C350_C01)
VlnPlot(C350_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C350_C01 <- subset(C350_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C350_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C350_C01)
## Pre-process Seurat object (sctransform)
C350_C01 <- SCTransform(C350_C01)
C350_C01 <- RunPCA(C350_C01)
C350_C01 <- RunUMAP(C350_C01, dims = 1:20)
####Doublet Finder
##C350_C01
## pK Identification (no ground-truth)
sweep.res.list_C350_C01 <- paramSweep_v3(C350_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C350_C01 <- summarizeSweep(sweep.res.list_C350_C01, GT = FALSE)
bcmvn_C350_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C350_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C350_C01 <- doubletFinder_v3(C350_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C350_C01 <- doubletFinder_v3(C350_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C350_C01@meta.data[,"DF_hi.lo"] <- C350_C01@meta.data$DF.classifications_0.25_0.09_82
C350_C01@meta.data$DF_hi.lo[which(C350_C01@meta.data$DF_hi.lo == "Doublet" & C350_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C350_C01@meta.data$DF_hi.lo[which(C350_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#350_C01
dim (C350_C01)
Idents(C350_C01) <- "DF_hi.lo"
C350_C01 <- subset(C350_C01, idents = "Singlet")
dim (C350_C01)
head (C350_C01@meta.data)
tail (C350_C01@meta.data)
####Remove unwanted genes (MT-)
#C350_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C350_C01)
total.genes <- list(rownames(C350_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C350_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C350_C01 <- subset (C350_C01, features = total.genes_subset)
dim(C350_C01)

####Pre-processing ----
####Load Dataset
C364_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C364_C01/outs/filtered_feature_bc_matrix")
####Create Seurat Object
C364_C01 <- CreateSeuratObject(counts = C364_C01.data, project = "C364_C01", min.cells = 0, min.features = 0)
#### Store mt.percent	
C364_C01 <- PercentageFeatureSet(C364_C01, pattern = "^MT-", col.name = "percent.mt")
###SUBSET low quality cells
dim (C364_C01)
VlnPlot(C364_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
C364_C01 <- subset(C364_C01, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
VlnPlot(C364_C01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
dim (C364_C01)
## Pre-process Seurat object (sctransform)
C364_C01 <- SCTransform(C364_C01)
C364_C01 <- RunPCA(C364_C01)
C364_C01 <- RunUMAP(C364_C01, dims = 1:20)
####Doublet Finder
##C364_C01
## pK Identification (no ground-truth)
sweep.res.list_C364_C01 <- paramSweep_v3(C364_C01, PCs = 1:15, sct = TRUE)
sweep.stats_C364_C01 <- summarizeSweep(sweep.res.list_C364_C01, GT = FALSE)
bcmvn_C364_C01 <- find.pK(sweep.stats_C402_C01)
## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(C364_C01@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies 
C364_C01 <- doubletFinder_v3(C364_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
C364_C01 <- doubletFinder_v3(C364_C01, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_82", sct = TRUE)
## PLOT the results 
C364_C01@meta.data[,"DF_hi.lo"] <- C364_C01@meta.data$DF.classifications_0.25_0.09_82
C364_C01@meta.data$DF_hi.lo[which(C364_C01@meta.data$DF_hi.lo == "Doublet" & C364_C01@meta.data$DF.classifications_0.25_0.09_82 == "Singlet")] <- "Doublet_lo"
C364_C01@meta.data$DF_hi.lo[which(C364_C01@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
##Remove doublets
#364_C01
dim (C364_C01)
Idents(C364_C01) <- "DF_hi.lo"
C364_C01 <- subset(C364_C01, idents = "Singlet")
dim (C364_C01)
head (C364_C01@meta.data)
tail (C364_C01@meta.data)
####Remove unwanted genes (MT-)
#C364_C01
#Remove -MT genes
####Get the total genes expressed list from seurat object
dim (C364_C01)
total.genes <- list(rownames(C364_C01))
total.genes <- as.data.frame(do.call(cbind, total.genes))
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = C364_C01), value = TRUE)
####get the MT-genes rows in the matrix
Mt.genes <- which(toupper(rownames(S402_C01)) %in% Mt.genes)
Mt.genes #### check which rows your MT- genes are
total.genes_subset <- total.genes[-Mt.genes,]
C364_C01 <- subset (C364_C01, features = total.genes_subset)
dim(C364_C01)


#### Merge objects----
Brain.FRONTAL = merge (C114_C01, y=c(C135_C01, C140_C01, C142_C01, C361_C01, C363_C01, C379_C01,C377_C01,
                                     C381_C01, C323_C01, C255_C01, C210_C01, C269_C01, C99_C01, C181_C01, C277_C01,
                                     C226_C01, C240_C01, C262_C01, C291_C01, C327_C01, C342_C01, C350_C01, C364_C01))
head (BrainFRONTAL@meta.data)
tail (BrainFRONTAL@meta.data)

####BrainFRONTALch - add metadata ----
write.csv(BrainFRONTAL@meta.data, file = "BrainFRONTAL_metadata.csv")
Brain.meta <- read.csv("BrainFRONTAL_metadata.csv")
rownames(Brain.meta)=colnames(BrainFRONTAL)
BrainFRONTAL <- AddMetaData(object=BrainFRONTAL, metadata=Brain.meta)

head(BrainFRONTAL@meta.data)
tail(BrainFRONTAL@meta.data)
.
.
.
.
.
.
.

#####Integrate data for Batch correction
head(BrainFRONTAL)
dim(BrainFRONTAL)
Idents (BrainFRONTAL)
Idents (BrainFRONTAL) <-"Patient"
VlnPlot(BrainFRONTAL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.0, ncol = 3)

###Check values of nFeature_RNA (number of genes)----
min<-min(BrainFRONTAL@meta.data$nFeature_RNA)
m<-median(BrainFRONTAL@meta.data$nFeature_RNA)
max<-max(BrainFRONTAL@meta.data$nFeature_RNA)
sd<-sd(BrainFRONTAL@meta.data$nFeature_RNA)

print(paste("Feature stats:",min,m,max,sd));

dim (BrainFRONTAL)

###Check values of nCount_RNA (number of UMI)----
min1<-min(BrainFRONTAL@meta.data$nCount_RNA)
m1<-median(BrainFRONTAL.combined@meta.data$nCount_RNA)
max1<-max(BrainFRONTAL@meta.data$nCount_RNA)
sd1<-sd(BrainFRONTAL@meta.data$nCount_RNA)

####nCount stats
Count95 <- quantile(BrainFRONTAL@meta.data$nCount_RNA, 0.95) # calculate value in the 95rd percentile

####Integrate to Perform Batch Correction----
FRONTAL.list <- SplitObject(BrainFRONTAL, split.by = "orig.ident")

#Normalize Data
FRONTAL.list <- lapply(X = FRONTAL.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = FRONTAL.list, nfeatures = 3000)
FRONTAL.list <- PrepSCTIntegration(object.list = FRONTAL.list, anchor.features = features)
Brain.anchors <- FindIntegrationAnchors(object.list = FRONTAL.list, normalization.method = "SCT", 
                                        anchor.features = features)
BrainFRONTAL.combined <- IntegrateData(anchorset = Brain.anchors, normalization.method = "SCT")

DefaultAssay(BrainFRONTAL.combined)
DefaultAssay(BrainFRONTAL.combined) <- "integrated"

##NOTE: To speed up this process for large datasets and maintain consistency as you collect new data, 
###can set one high-quality sample as a reference

#PCA Dimensionality Reduction
BrainFRONTAL.combined <- RunPCA(BrainFRONTAL.combined, features=VariableFeatures(object=BrainFRONTAL.combined))
ElbowPlot(BrainFRONTAL.combined)

#Visualization
BrainFRONTAL.combined <- RunUMAP(BrainFRONTAL.combined, reduction = "pca", dims = 1:30)
BrainFRONTAL.combined <- FindNeighbors(BrainFRONTAL.combined, dims = 1:30)
DimPlot(BrainFRONTAL.combined, group.by = "Heatmap.celltype", label = T)+NoLegend()
BrainFRONTAL.combined <- FindClusters(BrainFRONTAL.combined, resolution=6, verbose = FALSE,algorithm=1) 

BrainFRONTAL.combined@meta.data <-BrainFRONTAL.combined@meta.data [, -which(colnames(BrainFRONTAL.combined@meta.data) %in% 'integrated_snn_res.5')]

DimPlot (BrainFRONTAL.combined, label = T)
DimPlot (BrainFRONTAL.combined, label = F, split.by = "Diagnosis")+NoLegend()

BrainFRONTAL.combined <- NormalizeData(BrainFRONTAL.combined, normalization.method = "LogNormalize", scale.factor = 10000)

Thal <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$Thal.Phase)
BRAAK <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$BRAAK)
NIA.AA <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$NIA.AA)
AAD <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$AAD)
PMT <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$Post.mortem.time)
CERAD <- table(BrainFRONTAL.combined@active.ident, BrainFRONTAL.combined$CERAD)

##############################################x
##AddModule Score for cell type annotation----
##############################################x

AddModuleScore()

###Astrocytes
ASTSignature <- read.csv("cell_type_module_AST.csv", row.names = 1)
ASTSignature <- ASTSignature$Ast
ASTSignature <- list (ASTSignature)

###Oligodendrocytes
OLISignature <- read.csv("cell_type_module_OLI1.csv", row.names = 1)
OLISignature <- OLISignature$Oli
OLISignature <- list (OLISignature)

###MICROGLIA
MICSignature <- read.csv("cell_type_module_MIC.csv", row.names = 1)
MICSignature <- MICSignature$Mic
MICSignature <- list (MICSignature)

###OPCs
OPCSignature <- read.csv("cell_type_module_OPC.csv", row.names = 1)
OPCSignature <- OPCSignature$Opc
OPCSignature <- list (OPCSignature)

###ENDOTHELIAL
ENDOTSignature <- read.csv("cell_type_module_ENDOT_sel.csv", row.names = 1)
ENDOTSignature <- ENDOTSignature$Endot
ENDOTSignature <- list (ENDOTSignature)

ENDOTSignature1 <- list(c("EBF1", "ABCG2", "FLI1", "LEF1", "EMCN", 
                     "IFI27", "ADGRL4", "FLT1", "CLDN5"))

####PERICYTES
PERSignature <- read.csv("cell_type_module_PER.csv", row.names = 1)
PERSignature <- PERSignature$Per
PERSignature <- list (PERSignature)

###NEURONAL
NEUSignature <- read.csv("cell_type_module_NEU.csv", row.names = 1)
NEUSignature <- NEUSignature$Neu
NEUSignature <- list (NEUSignature)

DefaultAssay(BrainFRONTAL.combined) <- "RNA"

BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                     features = ASTSignature,
                                     ctrl = 5,
                                     name = 'AST')
BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = OLISignature,
                                        ctrl = 5,
                                        name = 'OLI')
BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = MICSignature,
                                        ctrl = 5,
                                        name = 'MIC')
BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = OPCSignature,
                                        ctrl = 5,
                                        name = 'OPC')
BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = PERSignature,
                                        ctrl = 5,
                                        name = 'PERI')
BrainFRONTAL.noMT <- AddModuleScore(object = BrainFRONTAL.noMT,
                                        features = ENDOTSignature,
                                        ctrl = 5,
                                        name = 'ENDOT')
BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = NEUSignature,
                                        ctrl = 5,
                                        name = 'NEUR')

BrainFRONTAL.combined <- AddModuleScore(object = BrainFRONTAL.combined,
                                        features = Excsignature,
                                        ctrl = 5,
                                        name = 'EXC')
DefaultAssay (BrainFRONTAL.combined) <- "RNA"

head (BrainFRONTAL.combined@meta.data)
Idents (BrainFRONTAL.combined) <- "integrated_snn_res.6"


VlnPlot(BrainFRONTAL.noMT, features = c("AST1", "OLI1", "MIC1", "OPC1", "PERI1", 
                                            "ENDOT1", "NEUR1", "Exc1", "Inh1"), 
        pt.size = 0, stack = T, flip = T, sort=T)+NoLegend()


###Ast cells----
####Astrocytes
Idents (BrainFRONTAL.noMT) <- "integrated_snn_res.6"
DimPlot (BrainFRONTAL.noMT)

Ast.cells <- subset (BrainFRONTAL.noMT, idents = c("15","16","21","39","44","53","55","66","70"))
DimPlot (Ast.cells)

DefaultAssay (Ast.cells) <- "integrated"

#Visualization
Ast.cells <- RunUMAP(Ast.cells, reduction = "pca", dims = 1:30)
Ast.cells <- FindNeighbors(Ast.cells, dims = 1:30)
Ast.cells <- FindClusters(Ast.cells, resolution=0.3, verbose = FALSE,algorithm=1) 

DimPlot(Ast.cells, label = T)
DefaultAssay (Ast.cells)
DefaultAssay (Ast.cells) <- "RNA"

Ast.cells_Diagnosis <- table (Ast.cells@active.ident, Ast.cells$Diagnosis)
write.csv(Ast.cells_Diagnosis, file = "Ast.cells_Diagnosis.csv")

Ast.cells_Patient <- table (Ast.cells@active.ident, Ast.cells$Patient)
write.csv(Ast.cells_Patient, file = "Ast.cells_Patient.csv")

Cluster.markers_Ast <- FindAllMarkers(object = Ast.cells, logfc.threshold = 0.25)
write.csv(Cluster.markers_Ast, file = "Cluster.markers_Ast.csv")

Idents (Ast.cells) <- "Diagnosis"
DefaultAssay(Ast.cells) <- "RNA"
dim (Ast.cells)

Ast_E280Avscontrol <- FindMarkers(Ast.cells, ident.1 = "E280A", 
                                   ident.2 = "Control", logfc.threshold = 0.0, min.pct = 0.25,
                                  max.cells.per.ident = 1005, random.seed = 50)
write.csv(Ast_E280Avscontrol=, file = "Ast_E280Avscontrol.csv")

Ast_SporadicvsControl <- FindMarkers(Ast.cells, ident.1 = "Sporadic", 
                                  ident.2 = "Control", logfc.threshold = 0.0, min.pct = 0.25,
                                  max.cells.per.ident = 1005, random.seed = 50)
write.csv(Ast_SporadicvsControl, file = "Ast_SporadicvsControl.csv")


###Mic cells----
####Microglia
head (BrainFRONTAL.noMT)
Idents (BrainFRONTAL.noMT) <- "integrated_snn_res.6"
Mic.cells <- subset (BrainFRONTAL.noMT, idents = c("20","23","30","88"))

DimPlot (BrainFRONTAL.noMT, label = T)+NoLegend()
DimPlot (Mic.cells, label = T)+NoLegend()
DimPlot (Mic.cells, label =T, split.by= "Diagnosis")+NoLegend()

DefaultAssay (Mic.cells) <- "integrated"

Mic_Patient <- table (Mic.cells@active.ident, Mic.cells@meta.data$Patient)
write.csv(Mic_Patient, file = " Mic.patient.csv")

#Visualization
Mic.cells <- RunUMAP(Mic.cells, reduction = "pca", dims = 1:30)
Mic.cells <- FindNeighbors(Mic.cells, dims = 1:30)

Mic.cells <- FindClusters(Mic.cells, resolution=0.2, verbose = FALSE,algorithm=1) 

DimPlot(Mic.cells, label = T)

DimPlot(Mic.cells, label = T, split.by="integrated_snn_res.0.2")+NoLegend()

Idents (Mic.cells) <- "integrated_snn_res.0.2"

cluster.id <- c("Mic-0", "Mic-1", "Mic-2")


names(cluster.id) <- levels(Mic.cells)
Mic.cells <- RenameIdents(Mic.cells, cluster.id)

Mic.cells[["subpop"]] <- Idents(object = Mic.cells)

DimPlot (Mic.cells, split.by = "Diagnosis")

DefaultAssay (Mic.cells)
DefaultAssay (Mic.cells) <- "RNA"

table (Mic.cells@active.ident, Mic.cells@meta.data$Diagnosis)

head (Mic.cells.integrated@meta.data)
DefaultAssay(Mic.cells)

####Cluster markers
Cluster.markers_Mic <- FindAllMarkers(object = Mic.cells, logfc.threshold = 0.25)
write.csv(Cluster.markers_Ast, file = "Cluster.markers_Ast.csv")


Mic_E280Avscontrol <- FindMarkers(Mic.cells, ident.1 = "E280A", 
                                  ident.2 = "Control", logfc.threshold = 0.0, min.pct = 0.25,
                                  max.cells.per.ident = 643, random.seed = 50)
write.csv(Mic_E280Avscontrol, file = "Mic_E280Avscontrol.csv")

Mic_Sporadicvscontrol <- FindMarkers(Mic.cells, ident.1 = "Sporadic", 
                                     ident.2 = "Control", logfc.threshold = 0.0, min.pct = 0.25)
write.csv(Mic_Sporadicvscontrol, file = "Mic_Sporadicvscontrol.csv")


###Oli cells----
####Oligodendrocytes
Oli.cells <- subset (BrainFRONTAL.noMT, idents = c("0","1","2","3","5","6", "8", "10","11",
                                                       "19","24","25","26","37","52","72","83"))
table (Oli.cells@active.ident)

dim(Oli.cells)
DimPlot (BrainFRONTAL.noMT, label = T)+NoLegend()
DimPlot (Oli.cells, label = T)+NoLegend()
DimPlot (Oli.cells, label =T, split.by= "Diagnosis")+NoLegend()

DefaultAssay (Oli.cells) <- "integrated"

#Visualization
Oli.cells <- RunUMAP(Oli.cells, reduction = "pca", dims = 1:30)
Oli.cells <- FindNeighbors(Oli.cells, dims = 1:30)
Oli.cells <- FindClusters(Oli.cells, resolution=0.4, verbose = FALSE,algorithm=1) 
DefaultAssay(Oli.cells)
DefaultAssay(Oli.cells) <- "RNA"

Idents (Oli.cells) <- "integrated_snn_res.0.2"
DimPlot (Oli.cells)

id <- c("Oli-0", "Oli-1", "Oli-2", "Oli-3", "Oli-4")

names(id) <- levels(Oli.cells)
Oli.cells <- RenameIdents(Oli.cells, id)
DimPlot (Oli.cells, label=T, repel = T)

Oli.cells[["subpop"]] <- Idents(object = Oli.cells)

table (Oli.cells@active.ident)

Oli.markers <- FindAllMarkers(Oli.cells, logfc.threshold = 0.0, min.pct=0.25)
write.csv(Oli.markers, file = "Oli_markers.csv")

write.csv(Oli0.markers, file = "Oli0_markers.csv")

DimPlot(Oli.cells, label = T)+NoLegend()
Oli_E280Avscontrol <- FindMarkers(Oli.cells, ident.1 = "E280A", 
                                  ident.2 = "Control", logfc.threshold = 0.25, min.pct = 0.25,
                                  max.cells.per.ident = 4430)
write.csv(Oli_E280Avscontrol, file = "Oli_E280Avscontrol.csv")

Oli_Sporadicvscontrol <- FindMarkers(Oli.cells.integrated, ident.1 = "Sporadic", 
                                     ident.2 = "Control", logfc.threshold = 0, min.pct = 0.25)
write.csv(Oli_Sporadicvscontrol, file = "Oli_Sporadicvscontrol.csv")

EnhancedVolcano(Oli_E280Avscontrol, x = 'avg_log2FC', y = 'p_val_adj', 
                pCutoff = 5e-2, FCcutoff = 0.25, lab = rownames(Oli_E280Avscontrol),
                title = "Oligodendrocytes", subtitle = "E280A vs Control", titleLabSize = 14, 
                subtitleLabSize = 12, captionLabSize = 12, axisLabSize = 12, legendLabSize = 12,
                xlim = c(-1,1))

EnhancedVolcano(Oli_Sporadicvscontrol, x = 'avg_log2FC', y = 'p_val_adj', 
                pCutoff = 5e-2, FCcutoff = 0.25, lab = rownames(Oli_Sporadicvscontrol),
                title = "Oligodendrocytes", subtitle = "Sporadic vs Control", titleLabSize = 14, 
                subtitleLabSize = 12, captionLabSize = 12, axisLabSize = 12, legendLabSize = 12,
                xlim = c(-1,1))


###Neu cells----
####Neuronal Cells
Neu.cells <- subset (BrainFRONTAL.complete, idents = c("0","1","2","3","5","6", "8", "10","11",
                                                   "19","24","25","26","37","52","72","83"))
table (Neu.cells@active.ident)

dim(Neu.cells)
DefaultAssay (Neu.cells) <- "integrated"

#Visualization
Neu.cells <- RunUMAP(Neu.cells, reduction = "pca", dims = 1:30)
Neu.cells <- FindNeighbors(Neu.cells, dims = 1:30)
Neu.cells <- FindClusters(Neu.cells, resolution=0.4, verbose = FALSE,algorithm=1) 
DefaultAssay(Neu.cells)
DefaultAssay(Neu.cells) <- "RNA"

Idents (Neu.cells) <- "integrated_snn_res.0.2"
DimPlot (Neu.cells)

Neu.markers <- FindAllMarkers(Neu.cells, logfc.threshold = 0.25, min.pct=0.25)
write.csv(Neu.markers, file = "Neu_markers.csv")


Exc_E280Avscontrol <- FindMarkers(Oli.cells, ident.1 = "E280A", 
                                  ident.2 = "Control", logfc.threshold = 0.25, min.pct = 0.25,
                                  max.cells.per.ident = 4430)
write.csv(Oli_E280Avscontrol, file = "Oli_E280Avscontrol.csv")

Oli_Sporadicvscontrol <- FindMarkers(Oli.cells.integrated, ident.1 = "Sporadic", 
                                     ident.2 = "Control", logfc.threshold = 0, min.pct = 0.25)
write.csv(Oli_Sporadicvscontrol, file = "Oli_Sporadicvscontrol.csv")

EnhancedVolcano(Oli_E280Avscontrol, x = 'avg_log2FC', y = 'p_val_adj', 
                pCutoff = 5e-2, FCcutoff = 0.25, lab = rownames(Oli_E280Avscontrol),
                title = "Oligodendrocytes", subtitle = "E280A vs Control", titleLabSize = 14, 
                subtitleLabSize = 12, captionLabSize = 12, axisLabSize = 12, legendLabSize = 12,
                xlim = c(-1,1))

EnhancedVolcano(Oli_Sporadicvscontrol, x = 'avg_log2FC', y = 'p_val_adj', 
                pCutoff = 5e-2, FCcutoff = 0.25, lab = rownames(Oli_Sporadicvscontrol),
                title = "Oligodendrocytes", subtitle = "Sporadic vs Control", titleLabSize = 14, 
                subtitleLabSize = 12, captionLabSize = 12, axisLabSize = 12, legendLabSize = 12,
                xlim = c(-1,1))

