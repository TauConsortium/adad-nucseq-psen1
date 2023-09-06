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
start_time <- Sys.time()

## Packages
library(lme4)
library(doParallel)
library(foreach)

################################################################################
## Functions
run_nbGLMM_for_gene <- function(gene_name, data_frame) {
    model <- glmer.nb(formula=paste0("`", gene_name, "` ~ Genotype + cdr + (1|Patient)"), 
                      data=data_frame)

    coefficients <- summary(model)$coefficients
    residuals <- residuals(model)
    aic <- AIC(model)

    # Return a list with the relevant results for the gene
    return(list(gene_name = gene_name,
                coefficients = coefficients,
                residuals = residuals,
                aic = aic))
}

################################################################################
## Inputs
# Define the number of cores to use for parallel processing
num_cores <- 10 

# outdirs
oDir <- "/home/eger/projects/Brain_snRNAseq/CCh/GLMMnb_E280A_6_Hets/"
cDir <- paste0(oDir, "coeff_results/")
rDir <- paste0(oDir, "residuals/")
aDir <- paste0(oDir, "aic/")

# choose cell type
ct <- "Ast"
ct_name <- sub(" ", "_", ct)

# load the .rds file
df <- readRDS(paste0(oDir, ct_name, "_df.RDS"))

################################################################################
## Run in parallel
# Initialize the parallel backend using doParallel
registerDoParallel(cores = num_cores)

# Get the gene names
gene_names <- colnames(df)[4:ncol(df)]

# Run the analysis in parallel for each gene using the same data frame
results_list <- foreach(gene_name = gene_names) %dopar% {
  run_nbGLMM_for_gene(gene_name, data_frame = df)
}

# Stop the parallel backend
stopImplicitCluster()

################################################################################
## Format and save results
# Extract all the coefficients and put them in single dataframe
coeff_df <- data.frame(Gene = gene_names, 
                        Estimate = sapply(results_list, function(x) x$coeff[2, 1]),
                        StdError = sapply(results_list, function(x) x$coeff[2, 2]),
                        z_value = sapply(results_list, function(x) x$coeff[2, 3]),
                        p_value = sapply(results_list, function(x) x$coeff[2, 4]))

# change to numeric
coeff_df[, 2:5] <- sapply(coeff_df[, 2:5], as.numeric)

# do multiple testing correction
coeff_df$fdr_p_value <- p.adjust(coeff_df$p_value, 'fdr')

# Extract all the residuals and put them in single dataframe
# column names are the gene names, rows are the cells
residuals_df <- data.frame(sapply(results_list, function(x) x$residuals))
colnames(residuals_df) <- gene_names

# Extract all the AIC values and put them in single dataframe 
aic_df <- data.frame(Gene = gene_names, 
                     AIC = sapply(results_list, function(x) x$aic))

# save all the dfs
out <- paste0(cDir, "glmer.nb1_", ct_name, "_coeff_results_hets.tsv")
write.table(coeff_df, file=out, quote=F, sep='\t', row.names=F)

out <- paste0(rDir, "glmer.nb1_", ct_name, "_residuals_hets.tsv")
write.table(residuals_df, file=out, quote=F, sep='\t')

out <- paste0(aDir, "glmer.nb1_", ct_name, "_aic_hets.tsv")
write.table(aic_df, file=out, quote=F, sep='\t', row.names=F)

################################################################################
end_time <- Sys.time()
# print elapsed time
print(end_time - start_time)
################################################################################
