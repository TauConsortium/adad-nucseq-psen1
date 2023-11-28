################################################################################
## Pipeline developed by Sarah J. Eger

## Step 6: Filter signicant genes
#   - must be significant & have Estimates with the same direction in 5 sample and 4 sample analyses

## Save significant gene lists
#  - for fisher's exact test
#  - for metascape upload
################################################################################

## Inputs

# best model
model_name <- "GLMMnb_indiv_cdr"

# tissues
cts <- c("WM", "GM")

# main outdir
oDir <- "/home/eger/projects/Brain_Visium/v5/"

outdir <- paste0(oDir, model_name, "/")
cDir <- paste0(outdir, "coeff_results/")

################################################################################
for (ct in cts){
    # 5 sample results
    fname1 <- paste0(cDir, ct, "_", model_name, "_coeff_results.tsv")
    coeff_df_5 <- read.table(fname1, sep="\t", header=TRUE)
    df1 <- coeff_df_5 # save copy

    # 4 sample results
    fname2 <- paste0(cDir, ct, "_4_sample_", model_name, "_coeff_results.tsv")
    coeff_df_4 <- read.table(fname2, sep="\t", header=TRUE)

    # print number of genes tested
    print(paste0(ct, " 4 samples total genes tested: ", nrow(coeff_df_4)))
    print(paste0(ct, " 5 samples total genes tested: ", nrow(coeff_df_5)))

    # remove non-significant genes
    coeff_df_4 <- coeff_df_4[coeff_df_4$fdr_p_value < 0.05,]
    coeff_df_5 <- coeff_df_5[coeff_df_5$fdr_p_value < 0.05,]

    # get lists of significantly upregulated or downregulated genes
    up_sig_4_genes <- coeff_df_4[coeff_df_4$Estimate > 0, "Gene"]
    down_sig_4_genes <- coeff_df_4[coeff_df_4$Estimate < 0, "Gene"]
    up_sig_5_genes <- coeff_df_5[coeff_df_5$Estimate > 0, "Gene"]
    down_sig_5_genes <- coeff_df_5[coeff_df_5$Estimate < 0, "Gene"]

    # print the length of each list
    print(paste0(ct, " 4 samples up genes: ", length(up_sig_4_genes)))
    print(paste0(ct, " 4 samples down genes: ", length(down_sig_4_genes)))
    print(paste0(ct, " 5 samples up genes: ", length(up_sig_5_genes)))
    print(paste0(ct, " 5 samples down genes: ", length(down_sig_5_genes)))
    
    # all sig genes in both analyses
    total_sig_4_5 <- unique(c(up_sig_4_genes, down_sig_4_genes, up_sig_5_genes, down_sig_5_genes))
    # all sig up genes in both analyses
    up_sig_4_5 <- unique(c(up_sig_4_genes, up_sig_5_genes))
    # all sig down genes in both analyses
    down_sig_4_5 <- unique(c(down_sig_4_genes, down_sig_5_genes))

    # sig genes with same direction in both analyses
    sig_4_5_overlap <- c(intersect(up_sig_4_genes, up_sig_5_genes), 
                        intersect(down_sig_4_genes, down_sig_5_genes))
    
    # print number of same upregulated genes out of 5 sample upregulated list
    print(paste0(length(intersect(up_sig_4_genes, up_sig_5_genes)), " out of ", 
            length(up_sig_5_genes), " upregulated ", ct, " genes in 5 samples are also upregulated in 4 samples"))
    
    # print number of same downregulated genes out of 5 sample downregulated list
    print(paste0(length(intersect(down_sig_4_genes, down_sig_5_genes)), " out of ", 
            length(down_sig_5_genes), " downregulated ", ct, " genes in 5 samples are also downregulated in 4 samples"))

    # create column for whether gene is in sig_4_5_overlap
    df1$Signicant_in_4_sample <- ifelse(df1$Gene %in% sig_4_5_overlap, "Yes", "No")

    # sort with 'Signicant_in_4_sample' 
    df1 <- df1[order(df1$fdr_p_value),]
    df1 <- df1[order(df1$Signicant_in_4_sample, decreasing=TRUE),]
    fname3 <- paste0(cDir, ct, "_", model_name, "_coeff_results_overlap.tsv")
    write.table(df1, fname3, sep="\t", row.names=FALSE, quote=FALSE)

    # in coeff_df_5, make any genes not in sig_4_5_overlap have a fdr_p_value of 1
    coeff_df_5[!(coeff_df_5$Gene %in% sig_4_5_overlap), "fdr_p_value"] <- 1

    # make XIST have a fdr_p_value of 1 because it is a sex effect
    coeff_df_5[coeff_df_5$Gene == "XIST", "fdr_p_value"] <- 1

    # save the results with the new fdr_p_values
    fname4 <- paste0(cDir, ct, "_", model_name, "_coeff_results_filtered.tsv")
    write.table(coeff_df_5, fname4, sep="\t", row.names=FALSE, quote=FALSE)

    # write the gene lists to files, no quotes
    up_file <- paste0(oDir, model_name, "/coeff_results/", ct, "_", model_name, "_up_genes_filtered.txt")
    down_file <- paste0(oDir, model_name, "/coeff_results/", ct, "_", model_name, "_down_genes_filtered.txt")
    write.table(intersect(up_sig_4_genes, up_sig_5_genes), up_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(intersect(down_sig_4_genes, down_sig_5_genes), down_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

