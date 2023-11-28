################################################################################
## Pipeline developed by Sarah J. Eger

## Step 5: Run a negative binomial GLMM on each gene
# Input: data for single cell type
# Output: GLMM results, AIC, residuals for each gene

# 5a_glmer.nb_WM_GM_prep.r
#   - prepares a dataframe for each type of matter in parallel
#   - only needs to be run once

# 5b_glmer.nb_WM_GM_parallelized_LME4.r
#   - runs the analysis for every gene in parallel for a single matter type
#   - needs to be run for each matter type (WM, GM)

# 5c_glmer.nb_4_sample_LME4.r
#   - same as 5b.glmer.nb_WM_GM_parallelized_LME4.r but with 4 samples
#   - C363_c1 is removed because it many more counts than the other samples
#   - needs to be run for each matter type (WM, GM)

################################################################################
## Conda environment
# conda activate Rann
# conda install -c conda-forge r-dharma
# conda install -c r r-foreach
# conda install -c conda-forge r-doparallel

################################################################################
## Packages
library(Seurat)
library(doParallel)
library(foreach)

################################################################################
## Functions
prepare_ct_dataframe <- function(ct, seurat_obj) {
    # stratify by cell type
    ct_obj <- subset(seurat_obj, MATTER == ct)

    # check number of samples
    if (length(unique(ct_obj@meta.data$orig.ident) == 4)){
        ct_name <- paste0(ct, "_4_sample")
    } else {
        ct_name <- ct
    }

    # set MATTER as identity
    Idents(ct_obj) <- ct_obj@meta.data$MATTER

    # get genes expressed in > 10% of nuclei
    dp <- DotPlot(ct_obj, features=rownames(ct_obj@assays$Spatial@counts))
    df <- dp$data
    rm("dp") # free up mem
    genes_exp <- rownames(df[df$pct.exp > 10, ])

    print(paste0(ct, " genes expressed: ", length(genes_exp)))

    counts <- t(as.matrix(ct_obj@assays$Spatial@counts[genes_exp, ]))
    df_counts <- as.data.frame(counts)
    df_counts[] <- sapply(df_counts, as.numeric)

    # get CDR
    cdr <- scale(rowSums(counts))

    # metadata cols needed
    df_meta <- ct_obj@meta.data[, c("INDIV", "REGION", "DX")]
    df_meta[] <- sapply(df_meta, as.factor)
    df_meta$cdr <- as.numeric(cdr)

    # combine
    df_m <- cbind(df_meta, df_counts)

    # force CN to be the reference level
    df_m <- within(df_m, DX <- relevel(factor(DX), ref = "CN"))

    # save the dataframe as an R object
    saveRDS(df_m, paste0(oDir, ct_name, "_df.RDS"))
}

################################################################################
## Inputs

# Define the number of cores to use for parallel processing
num_cores <- 2

# indir
iDir <- "/home/eger/projects/Brain_Visium/v3/clustering/"

# outdirs
oDir <- "/home/eger/projects/Brain_Visium/v4/GLMMnb_1/"

# celltypes = white/gray matter for spatial transcriptomics
cts <- c("WM", "GM")

################################################################################
## Format seurat object ##

# 5 samples - 3 HPC, 2 FR
samples <- c("C140_c12", "C363_c12", "C364_c12", "C363_c1",
             "C364_c1")

# load processed objects - get WM/GM
for (sample in samples) {
    load(paste0(iDir, sample, "/", sample, "_QC_cluster.Rdata"))
}

objs <- c(C140_c12, C363_c12, C364_c12, C363_c1,
            C364_c1)

# make new objects
for (obj in objs){
    sample <- obj@meta.data$orig.ident[1]
    df <- obj@meta.data[, c("seurat_clusters", "matter")]
    colnames(df) <- c("seurat_clusters", "MATTER")

    # load un-processed objects
    assign("qc_obj", get(load(paste0(iDir, sample, "/", sample, "_QC.Rdata"))))
    
    # add WM/GM assignments
    qc_obj <- AddMetaData(qc_obj, metadata = df)

    # remove spots not assigned to WM/GM
    qc_obj <- subset(qc_obj, MATTER != "Unknown")

    # remove genes detected in fewer than 2 spots in the sample
    all_counts <- qc_obj@assays$Spatial@counts
    all_counts <- all_counts[rowSums(all_counts != 0) >= 2,]
    qc_obj <- qc_obj[rownames(all_counts), ]

    # re-assign obj to original obj name
    assign(sample, qc_obj)
}

# merge samples
seurat_obj <- merge(C140_c12, y = c(C363_c12, C364_c12, C363_c1,
                C364_c1),
                add.cell.ids = samples, 
                project = "ALL")

################################################################################
## Format metadata ##

df <- seurat_obj@meta.data

# Diagnosis
df$DX <- "CN"
df$DX[df$orig.ident %in% c("C363_c12", "C140_c12", "C363_c1")] <- "E280A"

# Individual
df$INDIV <- df$orig.ident
df$INDIV[df$orig.ident %in% c("C363_c1", "C363_c12")] <- "C363"
df$INDIV[df$orig.ident %in% c("C364_c1", "C364_c12")] <- "C364"
df$INDIV[df$orig.ident == "C226_c12"] <- "C226"
df$INDIV[df$orig.ident == "C140_c12"] <- "C140"

# Brain region
df$REGION <- "FR"
df$REGION[df$orig.ident %in% c("C363_c12", "C140_c12", "C226_c12", "C364_c12")] <- "HPC"

seurat_obj <- AddMetaData(object = seurat_obj,
                    metadata = df)

# save object
saveRDS(seurat_obj, paste0(oDir, "Merged_5_samples_WM_GM.RDS"))

################################################################################
## Run in parallel
# Initialize the parallel backend using doParallel
registerDoParallel(cores = num_cores)

foreach(ct = cts) %dopar% {
    prepare_ct_dataframe(ct, seurat_obj)
}

# Stop the parallel backend
stopImplicitCluster()

################################################################################
## make 4 sample dataframes

seurat_obj <- readRDS(paste0(oDir, "Merged_5_samples_WM_GM.RDS"))

# remove C363_c1
seurat_obj <- subset(seurat_obj, orig.ident != "C363_c1")

# Initialize the parallel backend using doParallel
registerDoParallel(cores = num_cores)

foreach(ct = cts) %dopar% {
    prepare_ct_dataframe(ct, seurat_obj)
}

# Stop the parallel backend
stopImplicitCluster()
################################################################################