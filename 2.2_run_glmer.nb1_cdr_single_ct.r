################################################################################
# conda activate Rann
# conda install -c conda-forge r-dharma
################################################################################

library(Seurat)
library(lme4)
library(DHARMa)

# indir
iDir <- "/home/Camila/Sarah/"

# outdir
oDir <- "/home/eger/projects/Brain_snRNAseq/CCh/GLMMnb_E280A_3/"
dir.create(file.path(oDir), showWarnings = FALSE)

cDir <- paste0(oDir, "coeff_results/")
dir.create(file.path(cDir), showWarnings = FALSE)

rDir <- paste0(oDir, "residuals/")
dir.create(file.path(rDir), showWarnings = FALSE)

aDir <- paste0(oDir, "aic/")
dir.create(file.path(aDir), showWarnings = FALSE)

################################################################################
# load object
seurat_obj <- readRDS(paste0(iDir, "Brain.ch.RDS"))

# add cell type column
seurat_obj[["Celltype"]] <- Idents(object = seurat_obj)

# select E280As only (n = 8)
df <- seurat_obj@meta.data

# addmargins(table(df$Diagnosis_Reg, df$Celltype))

# select E280As only (n = 8)
df$ANALYSIS <- "YES"
df$ANALYSIS[(df$Diagnosis_Reg == "Sporadic") | 
                (df$Diagnosis_Reg == "Control") |
                (df$Diagnosis_Reg == "E280A ch OL")] <- "NO"
df$ANALYSIS[df$orig.ident == "C363_C01"] <- "NO"

# add label col
df$Genotype <- NULL
df$Genotype[df$Diagnosis_Reg == "E280A"] <- "E280A_WT"
df$Genotype[df$Patient == "E280A ch FC"] <- "E280A_Ch"

seurat_obj <- AddMetaData(object = seurat_obj,
                            metadata = df)

seurat_obj <- subset(seurat_obj, ANALYSIS == "YES")

################################################################################
# Only celltypes with enough nuclei
#cts <- c("Exc Neur", "Oli", "GABAergic Neur", "Ast", "Mic", "Endot/Peri")

# ct <- "Exc Neur"
# ct <- "Oli"
# ct <- "GABAergic Neur"
ct <- "Ast"

ct_name <- sub(" ", "_", ct)

# stratify by cell type
ct_obj <- subset(seurat_obj, Celltype == ct)

# get genes expressed in > 10% of nuclei
dp <- DotPlot(ct_obj, features=rownames(ct_obj@assays$RNA@counts))
df <- dp$data
rm("dp") # free up mem
genes_exp <- rownames(df[df$pct.exp > 10, ])

print(length(genes_exp))

counts <- as.matrix(t(ct_obj@assays$RNA@counts[genes_exp, ]))
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

# runn glmer.nb
glmer.nb1 <- lapply(genes_exp, function(x){
    list(gene = x,
            model = glmer.nb(formula=paste0("`", x, "` ~ Genotype + cdr + (1|Patient)"),
                data = df_m)
        )
    })

# get coeff and resid
outs <- lapply(glmer.nb1, function(x){
        list(gene = x$gene,
                coeff = summary(x$model)$coefficients,
                resid = residuals(x$model),
                aic = AIC(x$model)
        )
    })

# get p value
res <- as.data.frame(do.call(rbind, lapply(outs,
                                    function(i) c(
                                                Gene = i$gene[1], 
                                                Estimate = i$coeff[2, 1],
                                                StdError = i$coeff[2, 2],
                                                z_value = i$coeff[2, 3],
                                                p_value = i$coeff[2, 4]
                                                )
                                )
                    ))

# change to numeric
res[, 2:5] <- sapply(res[, 2:5], as.numeric)

# do multiple testing correction
res$fdr_p_value <- p.adjust(res$p_value, 'fdr')

# set identity to genotype
Idents(ct_obj) <- "Genotype"

# get avg_logFC and % of nuclei expressing
FCs <- FoldChange(ct_obj,
                ident.1 = "E280A_Ch",
                ident.2 = "E280A_WT",
                features = genes_exp,
                assay = "RNA")
colnames(FCs) <- c("avg_log2FC", "Ch_pct.exp", "WT_pct.exp")

# add FC
res <- cbind(FCs, res)

# reformat
res <- res[, c("Gene", "avg_log2FC", "Ch_pct.exp", "WT_pct.exp", 
                "Estimate", "StdError", "z_value", "p_value", "fdr_p_value")]
res <- res[order(res$fdr_p_value, -abs(res$avg_log2FC)), ]

# get resid df
resid <- as.data.frame(do.call(cbind, lapply(outs,
                                    function(i) i$resid)))
colnames(resid) <- genes_exp
resid <- cbind(df_meta, resid) # add metadata

# get AIC df
aic <- as.data.frame(do.call(rbind, lapply(outs,
                                    function(i) c(Gene = i$gene,
                                                    AIC = i$aic))))


# save all the dfs
out <- paste0(cDir, "glmer.nb1_", ct_name, "_coeff_results.tsv")
write.table(res, file=out, quote=F, sep='\t', row.names=F)

out <- paste0(rDir, "glmer.nb1_", ct_name, "_residuals.tsv")
write.table(resid, file=out, quote=F, sep='\t', row.names=F)

out <- paste0(aDir, "glmer.nb1_", ct_name, "_aic.tsv")
write.table(aic, file=out, quote=F, sep='\t', row.names=F)


