################################################################################
## Pipeline developed by Sarah J. Eger

## Step 7: Fisher's Exact Test for Spatial DEGs
# check for significant overlap between DEG lists

################################################################################
library(readxl)
library(tidyverse)

################################################################################
# Functions
# reformat snRNAseq results to match the spatial results
reformat_snRNAseq_GLMM_results <- function(df) {
    # rename 'gene' to 'Gene'
    colnames(df)[colnames(df) == "gene"] <- "Gene"

    # add a column for the 'fdr_p_value'
    df$fdr_p_value <- p.adjust(df$p_value, method="fdr")

    # remove the pseudobulk columns
    df <- df[, c("Gene", "Estimate", "StdError", "z_value", "p_value", "fdr_p_value")]

    # return the reformatted dataframe
    return(df)
}

# takes two dataframes with the results from two DGE analyses and checks if there is a significant overlap in the significant genes using a fisher's exact test
# returns a dataframe with the results of the fisher's exact test
dge_fisher_test <- function(overlap_name, df1, df2) {
    # default threshold
    threshold <- "fdr_p_value < 0.05"
    
    # get the list & number of significant genes in each dataframe
    sig1 <- df1$Gene[df1$fdr_p_value < 0.05]
    n_sig1 <- sum(df1$fdr_p_value < 0.05)

    sig2 <- df2$Gene[df2$fdr_p_value < 0.05]
    n_sig2 <- sum(df2$fdr_p_value < 0.05)

    # if there are no sig genes to compare to WM/GM, use the uncorrected p-value
    if (n_sig2 == 0) {
        sig2 <- df2$Gene[df2$p_value < 0.05]
        n_sig2 <- sum(df2$p_value < 0.05)
        threshold <- "p_value < 0.05"
    }

    # get the number of genes in each dataframe
    n_genes1 <- nrow(df1)
    n_genes2 <- nrow(df2)

    # get the number of genes that are significant in both dataframes
    n_sig_both <- length(intersect(sig1, sig2))

    # get the number of genes that are significant in one dataframe but not the other
    n_sig1_only <- length(setdiff(sig1, sig2))
    n_sig2_only <- length(setdiff(sig2, sig1))

    # get the number of genes that are not significant in either dataframe
    n_sig_none <- n_genes1 + n_genes2 - n_sig1 - n_sig2 + n_sig_both

    # create a dataframe with the results
    results_df <- data.frame(comparison = overlap_name,
                            threshold = threshold,
                            n_sig1 = n_sig1,
                            n_sig2 = n_sig2,
                            n_sig_both = n_sig_both,
                            n_sig1_only = n_sig1_only,
                            n_sig2_only = n_sig2_only,
                            n_sig_none = n_sig_none)

    # run the fisher's exact test
    fisher_result <- fisher.test(matrix(c(n_sig_both, n_sig1_only, n_sig2_only, n_sig_none), nrow=2))

    results_df$odds_ratio <- fisher_result$estimate
    results_df$conf_int_low <- fisher_result$conf.int[1]
    results_df$conf_int_high <- fisher_result$conf.int[2]
    results_df$p_value <- fisher_result$p.value

    # return the results dataframe
    return(results_df)
}


################################################################################
# spatial results
Dir1 <- "/home/eger/projects/Brain_Visium/v5/GLMMnb_indiv_cdr/coeff_results/"

# comparison data
Dir2 <- "/home/eger/projects/Brain_Visium/v5/DGE_comparisons/"

################################################################################
## Load all lists of DEGs needed

# Read in the results from the spatial analysis as tibbles
df_wm <- as_tibble(read.csv(paste0(Dir1, "WM_GLMMnb_indiv_cdr_coeff_results_filtered.tsv"), sep="\t"))
df_gm <- as_tibble(read.csv(paste0(Dir1, "GM_GLMMnb_indiv_cdr_coeff_results_filtered.tsv"), sep="\t"))

# Read in the GLMM results from the nuc-seq analysis
snRNAseq_results <- paste0(Dir2, "snRNAseq_DGEs_E280AvsControl.xlsx")
df_exc <- reformat_snRNAseq_GLMM_results(read_excel(snRNAseq_results, sheet="Exc"))
df_inh <- reformat_snRNAseq_GLMM_results(read_excel(snRNAseq_results, sheet="Inh"))
df_ast <- reformat_snRNAseq_GLMM_results(read_excel(snRNAseq_results, sheet="Ast"))
df_oli <- reformat_snRNAseq_GLMM_results(read_excel(snRNAseq_results, sheet="Oli"))

################################################################################
# get the results of the fisher's exact test
# WM & GM
df_wm_gm <- dge_fisher_test("WM & GM", df_wm, df_gm)
print(df_wm_gm)

# WM & Oli
df_wm_oli <- dge_fisher_test("WM & Oli", df_wm, df_oli)
print(df_wm_oli)

# GM & Exc
df_gm_exc <- dge_fisher_test("GM & Exc", df_gm, df_exc)
print(df_gm_exc)

# GM & Inh
df_gm_inh <- dge_fisher_test("GM & Inh", df_gm, df_inh)
print(df_gm_inh)

# combine all the results into one dataframe
df_all <- rbind(df_wm_gm, df_wm_oli, df_gm_exc, df_gm_inh)
print(df_all)

# write the results to a file
write.table(df_all, paste0(Dir2, "fisher_tests.tsv"), sep="\t", row.names=FALSE)

################################################################################