################################################################################
## Conda environment
# conda activate Rann
# conda install -c conda-forge r-dharma
# conda install -c r r-foreach
# conda install -c conda-forge r-doparallel

################################################################################
## Purpose: Run a negative binomial GLMM on each gene
# Input: data for single cell type
# Output: GLMM results, AIC, residuals for each gene

# 1_glmer.nb_all_cts_prep.r 
#   - prepares a dataframe for each ct in parallel
#   - only needs to be run once

# 2_glmer.nb_single_ct_parallelized.r 
#   - runs the analysis for every gene in parallel for a single ct
#   - needs to be run for each ct

################################################################################
## Packages
library(Seurat)
library(doParallel)
library(foreach)

################################################################################
## Functions
prepare_ct_dataframe <- function(ct, seurat_obj) {
    ct_name <- sub(" ", "_", ct)

    # stratify by cell type
    ct_obj <- subset(seurat_obj, Celltype == ct)

    # get genes expressed in > 10% of nuclei
    dp <- DotPlot(ct_obj, features=rownames(ct_obj@assays$RNA@counts))
    df <- dp$data
    rm("dp") # free up mem
    genes_exp <- rownames(df[df$pct.exp > 10, ])

    print(paste0(ct, " genes expressed: ", length(genes_exp)))

    counts <- t(as.matrix(ct_obj@assays$RNA@counts[genes_exp, ]))
    df_counts <- as.data.frame(counts)
    df_counts[] <- sapply(df_counts, as.numeric)

    # get CDR
    cdr <- scale(rowSums(counts))

    # metadata cols needed
    df_meta <- ct_obj@meta.data[, c("Patient", "Genotype")]
    df_meta[] <- sapply(df_meta, as.factor)
    df_meta$cdr <- as.numeric(cdr)

    # combine
    df_m <- cbind(df_meta, df_counts)

    # force WT to be the reference level
    df_m <- within(df_m, Genotype <- relevel(factor(Genotype), ref = "E280A_WT"))

    # save the dataframe as an R object
    saveRDS(df_m, paste0(oDir, ct_name, "_df.RDS"))
}

################################################################################
## Inputs

# Define the number of cores to use for parallel processing
num_cores <- 4 

# indir
iDir <- "/home/Camila/Sarah/"

# main object
main_obj <- readRDS(paste0(iDir, "BrainALL200.RDS"))

# outdirs
oDir <- "/home/eger/projects/Brain_snRNAseq/CCh/GLMMnb_E280A_6_Hets/"
cDir <- paste0(oDir, "coeff_results/")
rDir <- paste0(oDir, "residuals/")
aDir <- paste0(oDir, "aic/")

# create outdirs
for (dir in c(oDir, cDir, rDir, aDir)) {
    dir.create(file.path(dir), showWarnings = FALSE)
}

# celltypes
cts <- c("Exc Neur", "Oli", "GABAergic Neur", "Ast")

################################################################################
## Formatting metadata
# add cell type column
main_obj[["Celltype"]] <- Idents(object = main_obj)

# select E280As WTs and Hets only
df <- main_obj@meta.data

table(df$Diagnosis_Reg, df$Celltype)

# select E280As only (n = 8)
df$ANALYSIS <- "YES"
df$ANALYSIS[(df$Diagnosis_Reg == "Sporadic") | 
                (df$Diagnosis_Reg == "Control") |
                (df$Diagnosis_Reg == "E280A ch OL") | 
                (df$Diagnosis_Reg == "E280A ch FC")] <- "NO"

# add label col
df$Genotype <- NULL
df$Genotype[df$Diagnosis_Reg == "E280A"] <- "E280A_WT"
df$Genotype[df$Diagnosis_Reg == "E280A_ch_HET"] <- "E280A_Ch_HET"

table(df$Diagnosis_Reg, df$Genotype)
table(df$Diagnosis_Reg, df$ANALYSIS)

# add the metadata
main_obj <- AddMetaData(object = main_obj,
                            metadata = df)

# subset for analysis
seurat_obj <- subset(main_obj, ANALYSIS == "YES")

################################################################################
## Run in parallel
# Initialize the parallel backend using doParallel
registerDoParallel(cores = num_cores)

foreach(ct = cts) %dopar% {
    prepare_ct_dataframe(ct, seurat_obj)
}

# Stop the parallel backend
stopImplicitCluster()
