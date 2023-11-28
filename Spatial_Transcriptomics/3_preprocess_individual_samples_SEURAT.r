################################################################################

# Step 3: preprocess and cluster each sample individually with Seurat
## Pipeline developed by Sarah J. Eger based on: 
# https://satijalab.org/seurat/articles/spatial_vignette.html


## Set bioconda settings
## https://bioconda.github.io/
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --set channel_priority flexible


## Create env & Install packages
# conda create -n py39 python=3.9
# conda activate py39
# conda install -c conda-forge r-patchwork
# conda install -c bioconda/label/gcc7 r-seurat
# conda install -c conda-forge r-hdf5r
################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(hdf5r)

# spaceranger output directory
sDir <- "/home/eger/projects/Brain_Visium/spaceranger_outputs/"

# main directory
hDir <- "/home/eger/projects/Brain_Visium/v3/"
dir.create(file.path(hDir), showWarnings = FALSE)

# sub directories
cDir <- "/home/eger/projects/Brain_Visium/v3/clustering/"
dir.create(file.path(cDir), showWarnings = FALSE)
mDir <- "/home/eger/projects/Brain_Visium/v3/QC/"
dir.create(file.path(mDir), showWarnings = FALSE)

# List of 6 samples
samples <- c("C140_c12", "C226_c12", "C363_c12", "C364_c12",
             "C363_c1", "C364_c1")

# Quality Control Metrics
feats <- c("nFeature_Spatial", "nCount_Spatial", "percent_mito", 
            "percent_ribo", "percent_hb")

# loop thru samples to create single seuratobj
for (sample in samples){
    # Load image
    image_Dir <- paste0(sDir, sample, "/outs/spatial/")
    image <- Read10X_Image(image.dir = image_Dir, 
                        image.name = "tissue_lowres_image.png", 
                        filter.matrix = TRUE) # include spots that have been determined to be over tissue

    # Load data - col = spots row = features
    data_Dir <- paste0(sDir, sample, "/outs/")
    space_data <- Load10X_Spatial(data.dir = data_Dir, 
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial", 
                                slice = sample,
                                filter.matrix = TRUE,
                                image = image,
                                use.names = TRUE)

    # change name of object to sample ID
    Idents(space_data) <- sample
    space_data$orig.ident <- sample

    ### Visualize Quality Control on Spatial Plots ###
    space_data[["percent_mito"]] <- PercentageFeatureSet(space_data, pattern = "^MT-")
    space_data[["percent_ribo"]] <- PercentageFeatureSet(space_data, pattern = "^RP[SL]")
    space_data[["percent_hb"]] <- PercentageFeatureSet(space_data, pattern = "^HB[^(P)]")

    # make a spatial plot per sample per feature
    for (feat in feats) {
        png(paste0(mDir, sample, "_", feat, "_spatial.png"))
        
        # same ranges for all plots across samples
        if (feat == "nFeature_Spatial") {
            plot1 <- SpatialFeaturePlot(space_data, features = feat, 
                                        min.cutoff = 0, max.cutoff = 7500)
            } else if (feat == "nCount_Spatial") {
            plot1 <- SpatialFeaturePlot(space_data, features = feat, 
                                        min.cutoff = 0, max.cutoff = 40000)                
            } else {
            plot1 <- SpatialFeaturePlot(space_data, features = feat, 
                                        min.cutoff = 0, max.cutoff = 80) 
        } 
        print(plot1)
        dev.off()
    }

    # make a spatial plot for post-QC comparison
    space_data_sp <- space_data[, (space_data$percent_mito < 30) & 
                                (space_data$nFeature_Spatial > 100) & 
                                (space_data$nCount_Spatial > 0)]

    png(paste0(mDir, sample, "_spot_QC_spatial.png"))
    plot1 <- SpatialFeaturePlot(space_data, features = "percent_mito", 
                                min.cutoff = 0, max.cutoff = 80) 
    plot2 <- SpatialFeaturePlot(space_data_sp, features = "percent_mito", 
                                min.cutoff = 0, max.cutoff = 80) 
    plot3 <- wrap_plots(plot1, plot2) + 
                        plot_annotation(title = sample,
                        theme = theme(plot.title = element_text(hjust = 0.5, vjust = 0.5)))
    print(plot3)
    dev.off()

    # make the sample name the object variable name
    assign(sample, space_data)
}

################################################################################
### Visualize Quality Control on Violin Plots ###
# merge all the objs (only for plotting purposes)
all_data <- merge(C140_c12, y = c(C226_c12, C363_c12, C364_c12,
                C363_c1, C364_c1), 
                add.cell.ids = samples, project = "HPC_and_FRONTAL")

for (feat in feats){
    # all samples in single plot
    png(paste0(mDir, feat, "_violin.png"))
    plot1 <- VlnPlot(all_data, features = feat, pt.size = 0.1) +
                NoLegend() + 
                theme(axis.title.x=element_blank())
    print(plot1)
    dev.off()
}

################################################################################
### Filter spots ###
pre_table <- table(all_data$orig.ident)

objs <- c(C140_c12, C226_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    # less than 50% mito genes total, more than 100 features & more than 0 counts
    obj <- obj[, (obj$percent_mito < 30) & 
                (obj$nFeature_Spatial > 100) & 
                (obj$nCount_Spatial > 0)]
    
    # re-assign obj to original obj name
    assign(sample, obj)
}

all_data <- merge(C140_c12, y = c(C226_c12, C363_c12, C364_c12,
                C363_c1, C364_c1), 
                add.cell.ids = samples, project = "HPC_and_FRONTAL")

# check how many spots were removed
post_table <- table(all_data$orig.ident)

print(pre_table)
print(post_table)

# must recreate objs otherwise they are the unfiltered objects
objs <- c(C140_c12, C226_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

################################################################################
### Filter genes ###

# 1st check what are the top expressed genes in whole dataset
all_counts <- all_data@assays$Spatial@counts
all_counts@x <- all_counts@x/rep.int(colSums(all_counts), diff(all_counts@p))
most_expressed <- order(Matrix::rowSums(all_counts), decreasing = T)[20:1]
boxplot(t(as.matrix(all_counts[most_expressed, ])), cex = 0.1, las = 1,
    xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# Filter Mitocondrial genes
for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    iDir <- paste0(cDir, sample, "/")
    dir.create(file.path(iDir), showWarnings = FALSE)

    # remove mito genes
    obj <- obj[!grepl("^MT-", rownames(obj)), ]
    
    # re-assign obj to original obj name
    assign(sample, obj)

    # save QCed objects
    fname <- paste0(iDir, sample, "_QC.Rdata")
    save(list = sample, file = fname)

    # save list of spots and features that passed QC
    count_mtx <- obj@assays$Spatial@counts
    genes <- tibble(genes = rownames(count_mtx))
    spots <- tibble(spots = colnames(count_mtx))

    fname <- paste0(iDir, sample, "_QC_genes.txt")
    write.table(genes, file = fname, row.names = FALSE, quote=FALSE, col.names=FALSE)

    fname <- paste0(iDir, sample, "_QC_spots.txt")
    write.table(spots, file = fname, row.names = FALSE, quote=FALSE, col.names=FALSE)

}

# must recreate objs otherwise they are the unfiltered objects
objs <- c(C140_c12, C226_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

################################################################################
### NORMALIZATION & CLUSTERING ###

for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    iDir <- paste0(cDir, sample, "/")

    # run SCTransform
    obj <- SCTransform(obj, 
                        assay = "Spatial", 
                        method = "poisson")

    # Dimensionality reduction and clustering
    obj <- RunPCA(obj, 
                assay = "SCT")

    # Settings
    dimensions <- 1:10 # from Chen et al. (2021)
    res <- 0.1

    obj <- FindNeighbors(obj, 
                        reduction = "pca", 
                        dims = dimensions)
    obj <- FindClusters(obj,
                        resolution = res)
    obj <- RunUMAP(obj, 
                    reduction = "pca", 
                    dims = dimensions)

    # Number of clusters
    n <- length(levels(obj$seurat_clusters))

    print(sample)
    print(n)

    # re-assign obj to original obj name
    assign(sample, obj)
}

# must recreate objs otherwise they are the unfiltered objects
objs <- c(C140_c12, C226_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

################################################################################

for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    iDir <- paste0(cDir, sample, "/")

    # find markers for every cluster compared to all remaining cells
    # report only the positive ones
    markers <- FindAllMarkers(obj, 
                            only.pos = TRUE, # avg_log2FC direction
                            min.pct =  0.25) # only test genes that are detected in a minimum fraction of spots


    # save cluster markers to csv
    cols <- c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")
    markers <- markers[, cols]

    fname <- paste0(iDir, sample, "_cluster_top_genes.csv")
    write.csv(markers, file = fname, row.names = FALSE)

    # WM & GM cluster assignments
    GM.clus <- markers[markers$gene == "SNAP25", "cluster"]
    WM.clus <- markers[markers$gene == "MBP", "cluster"]

    # this sample does not clearly divide
    if (sample == "C226_c12"){
        df <- obj@meta.data
        df$matter <- "Unknown"

        obj <- AddMetaData(obj, metadata = df)
    } else {
        df <- obj@meta.data
        df$matter <- ifelse(df$seurat_clusters == GM.clus, "GM", "Unknown")
        df[df$seurat_clusters == WM.clus, "matter"] <- "WM"

        obj <- AddMetaData(obj, metadata = df)
    }

    print(sample)
    print(table(df$matter, df$seurat_clusters))

    # re-assign obj to original obj name
    assign(sample, obj)

    # save objects
    fname <- paste0(iDir, sample, "_QC_cluster.Rdata")
    save(list = sample, file = fname)

}
################################################################################