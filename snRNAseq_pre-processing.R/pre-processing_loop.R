
library(dplyr)
library(Seurat)
library (Matrix)
library(ggplot2)
library(sctransform)
library (EnhancedVolcano)
library (DoubletFinder)
library (pheatmap)


### Set up photo Directory
version <- 1
iDir <- "/home/daviswestover/projects/snRNAseq"
# create sub directory for version
vDir <- paste0(iDir, 'images', version, "/")
dir.create(vDir)

#Read in data
C114_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C114_C01/outs/filtered_feature_bc_matrix")
C135_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C135_C01/outs/filtered_feature_bc_matrix")
C140_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C140_C01/outs/filtered_feature_bc_matrix")
C142_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C142_C01/outs/filtered_feature_bc_matrix")
C361_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C361_C01/outs/filtered_feature_bc_matrix")
C363_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C363_C01/outs/filtered_feature_bc_matrix")
C379_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C379_C01/outs/filtered_feature_bc_matrix")
C377_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C377_C01/outs/filtered_feature_bc_matrix")
C381_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C381_C01/outs/filtered_feature_bc_matrix")
C323_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C323_C01/outs/filtered_feature_bc_matrix")
C255_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C255_C01/outs/filtered_feature_bc_matrix")
C210_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C210_C01/outs/filtered_feature_bc_matrix")
C269_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C269_C01/outs/filtered_feature_bc_matrix")
C99_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C99_C01/outs/filtered_feature_bc_matrix")
C181_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C181_C01/outs/filtered_feature_bc_matrix")
C277_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C277_C01/outs/filtered_feature_bc_matrix")
C226_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C226_C01/outs/filtered_feature_bc_matrix")
C240_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C240_C01/outs/filtered_feature_bc_matrix")
C262_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C262_C01/outs/filtered_feature_bc_matrix")
C291_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C291_C01/outs/filtered_feature_bc_matrix")
C327_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C327_C01/outs/filtered_feature_bc_matrix")
C342_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C342_C01/outs/filtered_feature_bc_matrix")
C350_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C350_C01/outs/filtered_feature_bc_matrix")
C364_C01.data <- Read10X(data.dir = "/home/Camila/BrainData/Raw_Data/Cellranger_outs/run_countpremRNA_C364_C01/outs/filtered_feature_bc_matrix"


####Create Seurat Objects
C114_C01 <- CreateSeuratObject(counts = C114_C01.data, project = "C114_C01", min.cells = 0, min.features = 0)
C135_C01 <- CreateSeuratObject(counts = C135_C01.data, project = "C135_C01", min.cells = 0, min.features = 0)
C140_C01 <- CreateSeuratObject(counts = C140_C01.data, project = "C140_C01", min.cells = 0, min.features = 0)
C142_C01 <- CreateSeuratObject(counts = C142_C01.data, project = "C142_C01", min.cells = 0, min.features = 0)
C361_C01 <- CreateSeuratObject(counts = C361_C01.data, project = "C361_C01", min.cells = 0, min.features = 0)
C363_C01 <- CreateSeuratObject(counts = C363_C01.data, project = "C363_C01", min.cells = 0, min.features = 0)
C379_C01 <- CreateSeuratObject(counts = C379_C01.data, project = "C379_C01", min.cells = 0, min.features = 0)
C377_C01 <- CreateSeuratObject(counts = C377_C01.data, project = "C377_C01", min.cells = 0, min.features = 0)
C381_C01 <- CreateSeuratObject(counts = C381_C01.data, project = "C381_C01", min.cells = 0, min.features = 0)
C323_C01 <- CreateSeuratObject(counts = C323_C01.data, project = "C323_C01", min.cells = 0, min.features = 0)
C255_C01 <- CreateSeuratObject(counts = C255_C01.data, project = "C255_C01", min.cells = 0, min.features = 0)
C210_C01 <- CreateSeuratObject(counts = C210_C01.data, project = "C210_C01", min.cells = 0, min.features = 0)
C269_C01 <- CreateSeuratObject(counts = C269_C01.data, project = "C269_C01", min.cells = 0, min.features = 0)
C99_C01 <- CreateSeuratObject(counts = C99_C01.data, project = "C99_C01", min.cells = 0, min.features = 0)
C181_C01 <- CreateSeuratObject(counts = C181_C01.data, project = "C181_C01", min.cells = 0, min.features = 0)
C277_C01 <- CreateSeuratObject(counts = C277_C01.data, project = "C277_C01", min.cells = 0, min.features = 0)
C226_C01 <- CreateSeuratObject(counts = C226_C01.data, project = "C226_C01", min.cells = 0, min.features = 0)
C240_C01 <- CreateSeuratObject(counts = C240_C01.data, project = "C240_C01", min.cells = 0, min.features = 0)
C262_C01 <- CreateSeuratObject(counts = C262_C01.data, project = "C262_C01", min.cells = 0, min.features = 0)
C291_C01 <- CreateSeuratObject(counts = C291_C01.data, project = "C291_C01", min.cells = 0, min.features = 0)
C327_C01 <- CreateSeuratObject(counts = C327_C01.data, project = "C327_C01", min.cells = 0, min.features = 0)
C342_C01 <- CreateSeuratObject(counts = C342_C01.data, project = "C342_C01", min.cells = 0, min.features = 0)
C350_C01 <- CreateSeuratObject(counts = C350_C01.data, project = "C350_C01", min.cells = 0, min.features = 0)
C364_C01 <- CreateSeuratObject(counts = C364_C01.data, project = "C364_C01", min.cells = 0, min.features = 0)

#Make list for loops
samples <- c(C114_C01, C135_C01, C140_C01, C142_C01, C361_C01, C363_C01, C379_C01, C377_C01,
C381_C01, C323_C01, C255_C01, C210_C01, C269_C01, C99_C01, C181_C01, C277_C01,
C226_C01, C240_C01, C262_C01, C291_C01, C327_C01, C342_C01, C350_C01, C364_C01)

#Empty list for pre-processed data
proc_samples <- c()


for (i in mini_samples) {
sample <- i
###Store Sample name
name <- sample@meta.data$orig.ident[1]
#### Store mt.percent
sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")


###SUBSET low quality cells
pre_sub_plot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
#Visualize low quality cells before subset
png(paste0(vDir, name, "_pre_featurecounts.png"))
print(pre_sub_plot)
dev.off() 

sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt < 5)
post_sub_plot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)  
#Visualize low quality cells after subset
png(paste0(vDir, name, "_post_featurecounts.png"))
print(post_sub_plot)
dev.off()   


## Pre-process Seurat object (sctransform)
sample <- SCTransform(sample)
sample <- RunPCA(sample)
sample <- RunUMAP(sample, dims = 1:20)

####Doublet Finder----
## pK Identification (no ground-truth)
sweep.res.list_sample <- paramSweep_v3(sample, PCs = 1:15, sct = TRUE)
sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
bcmvn_sample <- find.pK(sweep.stats_sample)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          	
nExp_poi <- round(0.030*length(sample@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Run DoubletFinder with varying classification stringencies 
sample <- doubletFinder_v3(sample, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
pANN_col_name <- colnames(sample@meta.data)[length(colnames(sample@meta.data))-1] # get caculated pANN value from meta.data
sample <- doubletFinder_v3(sample, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANN_col_name, sct = TRUE)



## Doublet Identifying Column; must be collected after running doublet finder
df_class_col_name <- colnames(sample@meta.data)[length(colnames(sample@meta.data))-1]

#Marking Doublets
sample@meta.data[,"DF_hi.lo"] <- sample@meta.data[[df_class_col_name]]
sample@meta.data$DF_hi.lo[which(sample@meta.data$DF_hi.lo == "Doublet" & sample@meta.data[[df_class_col_name]] == "Singlet")]  <- "Doublet_lo"
sample@meta.data$DF_hi.lo[which(sample@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"


##Remove doublets----
Idents(sample) <- "DF_hi.lo"
sample <- subset(sample, idents = "Singlet")



####Remove unwanted genes (MT-)----
#Remove -MT genes
####Get the total genes expressed list from seurat object
total.genes <- list(rownames(sample))
total.genes <- as.data.frame(do.call(cbind, total.genes))
#### Get all MT genes
Mt.genes <- grep(pattern = "^MT-", x = rownames(x = sample), value = TRUE)
#### Include S402_C01 as well
Mt.genes <- c(Mt.genes, 'S402_C01)')
#Subset out MT genes
total.genes_subset <- filter(total.genes, !(V1 %in% Mt.genes))
sample <- subset (sample, features = total.genes_subset$V1)
proc_samples <- c(proc_samples, sample)
}

#### Merge objects----
BrainFRONTAL = merge (proc_samples[[1]], y=c(proc_samples[[2]], proc_samples[[3]],proc_samples[[4]],proc_samples[[5]], 
proc_samples[[6]], proc_samples[[7]], proc_samples[[8]], proc_samples[[9]], proc_samples[[10]], proc_samples[[11]], proc_samples[[12]],
 proc_samples[[13]], proc_samples[[14]], proc_samples[[15]], proc_samples[[16]], proc_samples[[17]], proc_samples[[18]], proc_samples[[19]], 
 proc_samples[[20]], proc_samples[[21]] ,proc_samples[[22]], proc_samples[[23]], proc_samples[[24]]))

## Save merged/pre-processed data
saveRDS(proc_samples, file = "merged_samples.rds")
