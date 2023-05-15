################################################################################
# explore gene sets of interest from other papers

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
dDir <- "/home/eger/projects/Brain_Visium/v3/DGE_e/"
gDir <- "/home/eger/projects/Brain_Visium/v3/gene_sets/"

################################################################################
# modulators tau (and AB for some) aggregation
# https://www.science.org/doi/full/10.1126/sciadv.1600947
promoters <- c("FKBP5", "HS3ST2", "MAPT", "CDC37", "GSK3A", "GSK3B")
protectors <- c("HSPA1L", "HSPA1B", "HSPA1A", "BAG3", "DNAJA1", "HSPB8")
tangles <- read.csv(paste0(gDir, "tangles_Fu_2019.txt"))$Tangles

# Otero-Garcia 2022
# genes that are only up or down in at least one cluster
# genes both up & down in multiple clusters were removed
df_nft <- read.csv(paste0(gDir, "NFT_dysregulated_all.txt"), sep="\t") # 203
up_NFTs <- df_nft$Up # 162
down_NFTs <- df_nft$Down[1:41] # 41

# all upregulated in excitatory neurons in more than 3 clusters from their paper
NFTs_3 <- read.csv(paste0(gDir, "NFT_upregulated__genes_ranked.txt"), sep="\t")$gene # 95

# Camila's neuron genes - logthres = 0.1, pval adj < 0.05
df_nuc <- read.csv(paste0(gDir, "Neuronal_E280A_DGE.csv")) # 743

################################################################################
# first filter within our data
## GM DGE & NFT genes ##
reg <- "GM"

# DGE: 5 samples
fname <- fname <- paste0(dDir, "DGE_all_", reg, "_E280A_v_CN_5_samples.csv")
df1 <- read.csv(fname) # 6837

# DGE: 4 samples
fname <- paste0("/home/eger/projects/Brain_Visium/v3/DGE_d/DGE_all_", 
                reg, "_E280A_v_CN_4_samples.csv")
df2 <- read.csv(fname) # 4993

# sig filteriung
df1 <- df1[df1$p_val_adj < 0.05, ] # 5932 are significant p < 0.05
df1 <- df1[df1$p_val_adj < 0.05 & abs(df1$avg_log2FC) > 0.1, ]

df2 <- df2[df2$p_val_adj < 0.05, ] # 2767 are significant p < 0.05
df2 <- df2[df2$p_val_adj < 0.05 & abs(df2$avg_log2FC) > 0.1, ]

# only genes sig in both
df3 <- df1[df1$Gene %in% df2$Gene, ]

# shared direction genes
up_genes <- intersect(df1$Gene[df1$avg_log2FC > 0], df2$Gene[df2$avg_log2FC > 0]) # 1313
down_genes <- intersect(df1$Gene[df1$avg_log2FC < 0], df2$Gene[df2$avg_log2FC < 0]) # 213

df3 <- df3[df3$Gene %in% c(up_genes, down_genes), ] # 1526

################################################################################
# compare with NFT dataset

# significant in both datasets
dim(df3[df3$Gene %in% c(up_NFTs, down_NFTs), ]) # 67/203 are significant DEGs in GM

# upregulated in both datasets
NFT_up_spatial <- df3[df3$Gene %in% up_NFTs & df3$avg_log2FC > 0, "Gene"] # 45
length(NFT_up_spatial)

# downregulated in both datasets
NFT_down_spatial <- df3[df3$Gene %in% down_NFTs & df3$avg_log2FC < 0, "Gene"] # 2
length(NFT_down_spatial)

# save concordant genes to .csv
fname <- paste0(gDir, "NFT_concordant_GM_DGEs.csv")
write.csv(df3[df3$Gene %in% c(NFT_up_spatial, NFT_down_spatial), ], row.names=F, file=fname)

# 95 genes - only 39 are significant in the spatial data
NFT_up_3_spatial <- df3[df3$Gene %in% NFTs_3 & df3$avg_log2FC > 0, "Gene"] # 26
length(NFT_up_3_spatial)

################################################################################
# vol plot with 47 genes
reg <- "GM"

# DGE: 5 samples
fname <- fname <- paste0(dDir, "DGE_all_", reg, "_E280A_v_CN_5_samples.csv")
df1 <- read.csv(fname) # 6837

# DGE: 4 samples
fname <- paste0("/home/eger/projects/Brain_Visium/v3/DGE_d/DGE_all_", 
                reg, "_E280A_v_CN_4_samples.csv")
df2 <- read.csv(fname) # 4993

# shared direction genes
up_genes <- intersect(df1$Gene[df1$avg_log2FC > 0], df2$Gene[df2$avg_log2FC > 0]) # 2519
down_genes <- intersect(df1$Gene[df1$avg_log2FC < 0], df2$Gene[df2$avg_log2FC < 0]) # 1862 

dfg <- df2[df2$Gene %in% c(up_genes, down_genes), ]

keyvals <- ifelse(
dfg$avg_log2FC < -0.25, 'steelblue2',
    ifelse(dfg$avg_log2FC > 0.25, 'firebrick3',
    'gray100'))
keyvals[is.na(keyvals)] <- 'gray100'
names(keyvals)[keyvals == 'firebrick3'] <- 'Upregulated'
names(keyvals)[keyvals == 'steelblue2'] <- 'Downregulated'

p7 <- EnhancedVolcano(dfg, x = 'avg_log2FC', y = 'p_val_adj',

                        # labels
                        lab = dfg$Gene,
                        title = "",
                        subtitle = "",
                        caption = "",
                        selectLab = c(NFT_up_spatial, NFT_down_spatial),
                        colCustom = keyvals,
                        labSize = 5.5,
                        pointSize = 2,
                        boxedLabels = TRUE,
                        
                        # connectors
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,

                        # grid
                        gridlines.minor = FALSE,
                        gridlines.major = FALSE,                          
                        
                        # Axes
                        xlim = c(-2.5, 2.5),  
                        ylim = c(0, 310),        
                        xlab = bquote('Average' ~log[2]~ 'fold change'),
                        ylab = bquote('-' ~log[10]~ 'adjusted p-value'),
                        axisLabSize = 14,
                        
                        # legend
                        # legendIconSize = 3,
                        # legendPosition = "right",
                        # legendLabSize = 16,
                        
                        # thresholds
                        FCcutoff = 0.25,
                        pCutoff = 0.05,
                        ) +
    theme(axis.title=element_text(size=18),
            axis.line=element_line(size=1.5)) +
    guides(color = guide_legend(override.aes = list(size=3, alpha = 0.7)))  



################################################################################
## Spatial GM & snRNAseq neurons ##

# log thres = 0.1
#df4 <- df3[abs(df3$avg_log2FC) > 0.1, ] # 1718 are significant and log thres

# only genes sig in both
df4 <- df3[df3$Gene %in% df_nuc$Gene, ] # 143

# shared direction genes
up_genes <- intersect(df4$Gene[df4$avg_log2FC > 0.1], df_nuc$Gene[df_nuc$avg_log2FC > 0.1]) # 116
down_genes <- intersect(df4$Gene[df4$avg_log2FC < -0.1], df_nuc$Gene[df_nuc$avg_log2FC < -0.1]) # 10

df4 <- df4[df4$Gene %in% c(up_genes, down_genes), ] # 126

df4[df4$Gene %in% c(NFT_up_spatial, NFT_down_spatial), ]

#          Gene avg_log2FC pct.1 pct.2        p_val    p_val_adj
# 701      CNR1  0.1223594 0.566 0.277 1.810793e-56 4.257717e-52
# 846  CACNA2D1  0.2974803 0.672 0.366 4.405282e-53 1.035814e-48
# 1250    RAB10  0.1793336 0.685 0.387 1.745856e-45 4.105030e-41
# 1566   GABRA3  0.1755162 0.501 0.257 7.344780e-41 1.726978e-36
# 1729    NPTX1  0.2724702 0.907 0.651 5.572091e-39 1.310166e-34
# 2970   CNKSR2  0.1365994 0.745 0.484 2.925607e-27 6.878980e-23
# 3418    KCNJ3  0.1068142 0.623 0.385 3.273885e-24 7.697887e-20




# # check if genes overexpressed in cells with NFTs are upregulated in E280A gray matter
# dfg <- df1[df1$Gene %in% NFTs, ]

# df3 <- read.csv(paste0(gDir, "NFT_upregulated__genes_ranked.txt"), sep="\t")
# colnames(df3) <- c("Gene", "n AT8+ clusters expressing")

# df3 <- merge(df3, dfg, by="Gene")
# dim(df3[df3$avg_log2FC > 0 & df3$p_val_adj < 0.05, ]) # 45/95

# fname <- paste0(gDir, "NFT_upregulated_genes_gray_matter_DGE.csv")
# write.csv(df3, file=fname, row.names=F)

################################################################################

# # samples
# samples <- c("C140_c12", "C363_c12", "C364_c12", "C363_c1",
#              "C364_c1")

# # load processed objects - get WM/GM
# for (sample in samples) {
#     load(paste0(cDir, sample, "/", sample, "_QC_cluster.Rdata"))
# }

# objs <- c(C140_c12, C363_c12, C364_c12, C363_c1,
#             C364_c1)

# # make new objects
# for (obj in objs){
#     sample <- obj@meta.data$orig.ident[1]
#     df <- obj@meta.data[, c("seurat_clusters", "matter")]

#     # load un-processed objects
#     assign("qc_obj", get(load(paste0(cDir, sample, "/", sample, "_QC.Rdata"))))
    
#     # add WM/GM assignments
#     qc_obj <- AddMetaData(qc_obj, metadata = df)

#     # remove genes detected in fewer than 2 spots in the sample
#     all_counts <- qc_obj@assays$Spatial@counts
#     all_counts <- all_counts[rowSums(all_counts != 0) >= 2,]
#     qc_obj <- qc_obj[rownames(all_counts), ]

#     # re-assign obj to original obj name
#     assign(sample, qc_obj)
# }

# SpatialFeaturePlot(C363_c12, features=promoters)

