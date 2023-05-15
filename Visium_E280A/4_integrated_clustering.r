################################################################################
# make UMAP with SCT data
# https://smorabit.github.io/hdWGCNA/articles/ST_basics.html
# https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette.html

# conda activate py39
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)

# sub directories
cDir <- "/home/eger/projects/Brain_Visium/v3/clustering/"
pDir <- "/home/eger/projects/Brain_Visium/v3/PCA/"
dir.create(file.path(pDir), showWarnings = FALSE)

# List of 5 samples
samples <- c("C140_c12",  "C363_c12", "C364_c12",
             "C363_c1", "C364_c1")

# load processed objects
for (sample in samples) {
    load(paste0(cDir, sample, "/", sample, "_QC_cluster.Rdata"))
}

objs <- c(C140_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

# load un-processed objects & add metadata
for (obj in objs){
    sample <- obj@meta.data$orig.ident[1]
    df <- obj@meta.data[, c("seurat_clusters", "matter")]

    assign("qc_obj", get(load(paste0(cDir, sample, "/", sample, "_QC.Rdata"))))
    qc_obj <- AddMetaData(qc_obj, metadata = df)

    # re-assign obj to original obj name
    assign(sample, qc_obj)

}
################################################################################
### INTEGRATION ###

# create a list of samples
sample_list <- list(C140_c12 = C140_c12,
                C363_c12 = C363_c12, C364_c12 = C364_c12,
                C363_c1 = C363_c1, C364_c1 = C364_c1)

# run SCT on all samples
sample_list <- lapply(sample_list, SCTransform, assay = "Spatial", method = "poisson")


# Need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB

# Select 5,000 features for integrating - 3k only has 5 CA1 genes
sample_features <- SelectIntegrationFeatures(sample_list, 
                                        nfeatures = 5000, 
                                        verbose = TRUE)

sample_list <- PrepSCTIntegration(object.list = sample_list, 
                            anchor.features = sample_features,
                            verbose = TRUE)

# Perform integration
sample_int.anchors <- FindIntegrationAnchors(object.list = sample_list, 
                                        normalization.method = "SCT",
                                        anchor.features = sample_features,
                                        verbose = TRUE)

sample_integrated <- IntegrateData(anchorset = sample_int.anchors, 
                                normalization.method = "SCT",
                                verbose = TRUE)

# Run dimensionality reduction and clustering
sample_integrated <- RunPCA(sample_integrated, verbose = TRUE)

# Check elbow plot
ElbowPlot(sample_integrated)

sample_integrated <- FindNeighbors(sample_integrated, dims = 1:10)

# clustering
res <- 0.1
sample_integrated <- FindClusters(sample_integrated, 
                            resolution = res,
                            verbose = FALSE)
sample_integrated <- RunUMAP(sample_integrated, reduction = "pca", dims = 1:10)

################################################################################
# add metadata
df <- sample_integrated@meta.data
df$DX <- "CN"
df$DX[df$orig.ident %in% c("C363_c12", "C140_c12", "C363_c1")] <- "E280A"

df$SUBREGION <- "NA"
df$SUBREGION[(df$REGION == "HPC") & (df$matter == "GM")] <- "HPC_GM"
df$SUBREGION[(df$REGION == "HPC") & (df$matter == "WM")] <- "HPC_WM"
df$SUBREGION[(df$REGION == "FR") & (df$matter == "GM")] <- "FR_GM"
df$SUBREGION[(df$REGION == "FR") & (df$matter == "WM")] <- "FR_WM"

df$matter <- ifelse(df$matter == "WM", "White Matter", df$matter)
df$matter[df$matter == "GM"] <- "Gray Matter"

df$REGION <- "Frontal Cortex"
df$REGION[df$orig.ident %in% c("C363_c12", "C140_c12", "C226_c12", "C364_c12")] <- "Hippocampus"

df$INDIV <- df$orig.ident
df$INDIV[df$orig.ident %in% c("C363_c1", "C363_c12")] <- "C363"
df$INDIV[df$orig.ident %in% c("C364_c1", "C364_c12")] <- "C364"
df$INDIV[df$orig.ident == "C226_c12"] <- "C226"
df$INDIV[df$orig.ident == "C140_c12"] <- "C140"

sample_integrated <- AddMetaData(object = sample_integrated,
                    metadata = df)

# save objects
save(list = c("sample_integrated"), file = paste0(cDir,"sample_integrated.Rdata"))

################################################################################
### PLOTS ###

Dplot1 <- DimPlot(sample_integrated, 
                reduction = "umap",
                group.by = "ident",
                )

Dplot2 <- DimPlot(sample_integrated, 
                reduction = "umap",
                group.by = "orig.ident",
                )

Dplot3 <- DimPlot(sample_integrated, 
                reduction = "umap",
                group.by = "matter",
                cols = c("#A6CEE3", "#B2DF8A", "#1F78B4")) +
                theme(legend.position = c(0.7, 0.2),
                legend.text=element_text(size=15))

Dplot4 <- DimPlot(sample_integrated, 
                reduction = "umap",
                group.by = "REGION",
                cols = c("#66C2A5", "#FC8D62")) +
                theme(legend.position = c(0.7, 0.2),
                legend.text=element_text(size=15))

# Idents(C363_c1) <- C363_c1@meta.data$matter
# SpatialDimPlot(C363_c1,
#             #label.size=10,
#             label.box = F) +
#             scale_fill_brewer(palette = "Paired") +
#             theme(legend.position = c(0.12, 0.05), #R/L, U/D
#             legend.text=element_text(size=15),
#             legend.title = element_blank())


# # check PC loadings and plots
# Lplot <- VizDimLoadings(sample_integrated, 
#                         dims = 1:5, 
#                         reduction = "pca",
#                         ncol = 5)
# print(Lplot)

# Dplot1 <- DimPlot(sample_integrated, 
#                 reduction = "umap",
#                 group.by = "REGION",
#                 cols = c("#66C2A5", "#FC8D62"))

# Dplot2 <- DimPlot(sample_integrated, 
#                 reduction = "pca",
#                 group.by = "matter",
#                 cols = c("#A6CEE3", "#1F78B4", "#B2DF8A"))

# Dplot3 <- DimPlot(ssample_integrated, 
#                 reduction = "pca",
#                 group.by = "DX",
#                 dims = c(1, 3))


# Dplot <- DimPlot(seurat_obj, 
#                 reduction = "pca",
#                 group.by = "SUBREGION")

# Dplot <- DimPlot(seurat_obj, 
#                 reduction = "pca",
#                 group.by = "INDIV")






