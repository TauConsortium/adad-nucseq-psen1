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
start_time <- Sys.time()

## Packages
library(lme4)
library(MASS)
library(doParallel)
library(foreach)

################################################################################
## Functions
run_model_for_gene <- function(gene_name, formula_str, data_frame) {
    print(paste0("Running model for ", gene_name, " with formula ", formula_str))

    model <- glmer.nb(formula=paste0("`", gene_name, "` ", formula_str), 
                        data=data_frame,
                        nAGQ=0)

    coefficients <- summary(model)$coefficients
    residuals <- residuals(model)
    aic <- AIC(model)

    # Return a list with the relevant results for the gene
    return(list(gene_name = gene_name,
                coefficients = coefficients,
                residuals = residuals,
                aic = aic))
}

compile_results_for_model <- function(results_list){
    # Extract all the coefficients and put them in single dataframe
    coeff_df <- data.frame(Gene = gene_names, 
                            Estimate = sapply(results_list, function(x) x$coefficients[2, 1]),
                            StdError = sapply(results_list, function(x) x$coefficients[2, 2]),
                            z_value = sapply(results_list, function(x) x$coefficients[2, 3]),
                            p_value = sapply(results_list, function(x) x$coefficients[2, 4]))

    # change to numeric
    coeff_df[, 2:5] <- sapply(coeff_df[, 2:5], as.numeric)

    # do multiple testing correction
    coeff_df$fdr_p_value <- p.adjust(coeff_df$p_value, 'fdr')

    # Extract all the residuals and put them in single dataframe
    # column names are the gene names, rows are the cells
    residuals_df <- data.frame(sapply(results_list, function(x) x$residuals))
    colnames(residuals_df) <- gene_names

    # Extract all the AICs and put them in single dataframe
    aic_df <- data.frame(Gene = gene_names, 
                        AIC = sapply(results_list, function(x) x$aic))

    # return the dataframes of the results
    return(list(coeff_df = coeff_df,
                residuals_df = residuals_df,
                aic_df = aic_df))
}

################################################################################
## Inputs
# Define the number of cores to use for parallel processing
num_cores <- 20

model_name <- "GLMMnb_indiv_cdr"
formula_str <- "~ DX + cdr + (1|INDIV)"

# choose cell type
ct <- "WM"
ct_name <- paste0(ct, "_4_sample")

# main outdir
oDir <- "/home/eger/projects/Brain_Visium/v5/"

outdir <- paste0(oDir, model_name, "/")
cDir <- paste0(outdir, "coeff_results/")
rDir <- paste0(outdir, "residuals/")
aDir <- paste0(outdir, "aic/")

# load the corresponding .rds file
df <- readRDS(paste0("/home/eger/projects/Brain_Visium/v4/GLMMnb_1/", ct_name, "_df.RDS"))

################################################################################
## Run in parallel
# Initialize the parallel backend using doParallel
registerDoParallel(cores = num_cores)

# Get the gene names
gene_names <- colnames(df)[5:ncol(df)]


# loop through the genes
model_results <- foreach(gene_name = gene_names) %dopar% {
    run_model_for_gene(gene_name, formula_str = formula_str, data_frame = df)
}

# compile the results
compiled_results <- compile_results_for_model(model_results)

# save the compiled results as 3 tsv files
write.table(compiled_results$coeff_df, paste0(cDir, ct_name, "_", model_name, "_coeff_results.tsv"), sep="\t", row.names=FALSE)
write.table(compiled_results$residuals_df, paste0(rDir, ct_name, "_", model_name, "_residuals.tsv"), sep="\t", row.names=FALSE)
write.table(compiled_results$aic_df, paste0(aDir, ct_name, "_", model_name, "_aic.tsv"), sep="\t", row.names=FALSE)

stopImplicitCluster()

################################################################################
end_time <- Sys.time()
# print elapsed time
print(end_time - start_time)
################################################################################
