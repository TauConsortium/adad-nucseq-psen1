################################################################################
# Plots for figure 1

# conda activate py39
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(EnhancedVolcano)

# sub directories
cDir <- "/home/eger/projects/Brain_Visium/v3/clustering/"
dDir <- "/home/eger/projects/Brain_Visium/v3/DGE_e/"
vDir <- "/home/eger/projects/Brain_Visium/v3/plots/"
dir.create(file.path(vDir), showWarnings = FALSE)

# List of 6 samples
samples <- c("C140_c12", "C226_c12", "C363_c12", "C364_c12",
             "C363_c1", "C364_c1")

# load objects
for (sample in samples) {
    load(paste0(cDir, sample, "/", sample, "_QC_cluster.Rdata"))
}

# integrated obj
load(paste0(cDir,"sample_integrated.Rdata"))

objs <- c(C140_c12, C226_c12, C363_c12, C364_c12, 
            C363_c1, C364_c1)

################################################################################
## PANEL A & B ##

# Spatial SCT - SNAP25 & MBP
for (gene in c("SNAP25", "MBP")){
    p1 <- SpatialFeaturePlot(C363_c1,
                            alpha =  c(0.3, 1), # min, max
                            image.alpha = 0.4,
                            crop = FALSE,
                            features = gene) + 
                            scale_fill_distiller(
                                type = "seq", 
                                palette = "BuPu", 
                                direction=1,
                                guide = guide_colourbar(direction = "vertical",
                                                    frame.colour = "black",
                                                    frame.linewidth = 0.7)) +
                            theme(legend.position=c(0.09, 0.13), # x, y
                                legend.text=element_text(size=15),
                                legend.title=element_blank(),
                                legend.background = element_rect(fill = NA))

    png(paste0(vDir, "Fig1a_", gene, ".png"))
    print(p1)
    dev.off()
}

################################################################################
## PANEL C, D, E ##

# change dx labels
df <- sample_integrated@meta.data
df$DX[df$DX == "CN"] <- "Control"
sample_integrated <- AddMetaData(sample_integrated, metadata = df)

# UMAPs for WM/GM
Idents(sample_integrated) <- as.factor(sample_integrated@meta.data$matter)

p2 <- DimPlot(sample_integrated, 
            reduction = "umap",
            cols = c("Gray Matter" = "#00C1AA", 
                    "White Matter" = "#FC61C9",
                    "Unknown" = "gray35")) +
            theme(legend.position=c(0.7, 0.93), # x, y
                legend.text=element_text(size=18),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                axis.text.x=element_text(size=14),
                axis.text.y=element_text(size=14),
                axis.line = element_line(size = 1),
                axis.ticks = element_line(size = 1)) +
            scale_x_continuous(limits = c(-10, 10)) +
            scale_y_continuous(limits = c(-5, 5),
                breaks=c(-3, 0, 3))

png(paste0(vDir, "Fig1b_matter.png"))
print(p2)
dev.off()

# UMAPs for HPC/FR
Idents(sample_integrated) <- as.factor(sample_integrated@meta.data$REGION)

p3 <- DimPlot(sample_integrated, 
            reduction = "umap",
            cols = c("Hippocampus" = "#E7861B", 
                    "Frontal Cortex" = "#AC88FF")) +
            theme(legend.position=c(0.7, 0.93), # x, y
                legend.text=element_text(size=18),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                axis.text.x=element_text(size=14),
                axis.text.y=element_text(size=14),
                axis.line = element_line(size = 1),
                axis.ticks = element_line(size = 1)) +
            scale_x_continuous(limits = c(-10, 10)) +
            scale_y_continuous(limits = c(-5, 5),
                breaks=c(-3, 0, 3))

png(paste0(vDir, "Fig1b_region.png"))
print(p3)
dev.off()

# UMAPs for DX
Idents(sample_integrated) <- as.factor(sample_integrated@meta.data$DX)

p4 <- DimPlot(sample_integrated, 
            reduction = "umap",
            cols = c("Control" = "#F8766D", 
                    "E280A" = "#00B81F")) +
            theme(legend.position=c(0.7, 0.93), # x, y
                legend.text=element_text(size=18),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                axis.text.x=element_text(size=14),
                axis.text.y=element_text(size=14),
                axis.line = element_line(size = 1),
                axis.ticks = element_line(size = 1)) +
            scale_x_continuous(limits = c(-10, 10)) +
            scale_y_continuous(limits = c(-5, 5),
                breaks=c(-3, 0, 3))

png(paste0(vDir, "Fig1b_dx.png"))
print(p4)
dev.off()


################################################################################
## PANEL F ##

# change matter labels
df <- C363_c1@meta.data

df$matter[df$matter == "GM"] <- "Gray Matter"
df$matter[df$matter == "WM"] <- "White Matter"

C363_c1 <- AddMetaData(C363_c1, metadata = df)

Idents(C363_c1) <- as.factor(C363_c1@meta.data$matter)

p5 <- SpatialDimPlot(C363_c1,
                alpha = c(0.5, 1), # min, max
                image.alpha = 0.4,
                crop = FALSE,
                cols = c("Gray Matter" = "#00C1AA", 
                            "White Matter" = "#FC61C9")) +
                theme(legend.position=c(0.2, 0.13), # x, y
                    legend.text=element_text(size=16),
                    legend.title=element_blank(),
                    legend.background = element_rect(fill = NA),
                    legend.key = element_blank()) +
                scale_size_manual() +
                guides(color = guide_legend(override.aes = list(
                    size = 10,
                    alpha = 1)))

png(paste0(vDir, "Fig1c_spatial.png"))
print(p5)
dev.off()                



    #             scale_fill_manual(values = c("#002F70", "#EDB4B5")) +
    #             guides(color = guide_legend(override.aes = list(
    #                 size = 10)))

    # geom_point(alpha = 0, size = 5) +
    # guides(color = guide_legend(override.aes = list(
    #     size = 10,
    #     alpha = 1)))

################################################################################
# volcano plots

# genes referenced in the text
txt_genes <- c("HSPA1A", "GJA1", "APOE", "CLU", "TXNIP")

################################################################################
reg <- "WM"

# DGE: 5 samples
fname <- fname <- paste0(dDir, "DGE_all_", reg, "_E280A_v_CN_5_samples.csv")
df1 <- read.csv(fname) # 

# DGE: 4 samples
fname <- paste0("/home/eger/projects/Brain_Visium/v3/DGE_d/DGE_all_", 
                reg, "_E280A_v_CN_4_samples.csv")
df2 <- read.csv(fname) # 

# shared direction genes
up_genes <- intersect(df1$Gene[df1$avg_log2FC > 0], df2$Gene[df2$avg_log2FC > 0]) # 
down_genes <- intersect(df1$Gene[df1$avg_log2FC < 0], df2$Gene[df2$avg_log2FC < 0]) # 

dfw <- df1[df1$Gene %in% c(up_genes, down_genes), ]

# genes to label
top_wm_genes <- dfw[order(-abs(dfw$avg_log2FC)), "Gene"][1:10]
top_wm_genes <- unique(c(top_wm_genes, txt_genes))

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
dfw$avg_log2FC < -0.25, 'steelblue2',
    ifelse(dfw$avg_log2FC > 0.25, 'firebrick3',
    'ghostwhite'))
keyvals[is.na(keyvals)] <- 'ghostwhite'
names(keyvals)[keyvals == 'firebrick3'] <- 'Upregulated'
names(keyvals)[keyvals == 'steelblue2'] <- 'Downregulated'

p6 <- EnhancedVolcano(dfw, x = 'avg_log2FC', y = 'p_val_adj',
                        # labels
                        lab = dfw$Gene,
                        title = "",
                        subtitle = "",
                        caption = "",
                        # boxedLabels = TRUE,
                        # selectLab = c("APOE", "AQP4"),   
                        selectLab = top_wm_genes,
                        colCustom = keyvals,
                        labSize = 5.5,
                        pointSize = 2,
                        
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
                        legendLabSize = 16,
                        
                        # thresholds
                        FCcutoff = 0.25,
                        pCutoff = 0.05,
                        ) +
    theme(legend.position=c(0.85, 0.2),
            axis.title=element_text(size=18),
            axis.line=element_line(size=1.5)) +
    guides(color = guide_legend(override.aes = list(size=3, alpha = 0.7))) 

png(paste0(vDir, "Fig1_WM_Volcano.png"), width = 600, height = 550)
print(p6)
dev.off()

################################################################################

reg <- "GM"

# DGE: 5 samples
fname <- fname <- paste0(dDir, "DGE_all_", reg, "_E280A_v_CN_5_samples.csv")
df1 <- read.csv(fname) # 

# DGE: 4 samples
fname <- paste0("/home/eger/projects/Brain_Visium/v3/DGE_d/DGE_all_", 
                reg, "_E280A_v_CN_4_samples.csv")
df2 <- read.csv(fname) # 

# shared direction genes
up_genes <- intersect(df1$Gene[df1$avg_log2FC > 0], df2$Gene[df2$avg_log2FC > 0]) # 
down_genes <- intersect(df1$Gene[df1$avg_log2FC < 0], df2$Gene[df2$avg_log2FC < 0]) # 

dfg <- df1[df1$Gene %in% c(up_genes, down_genes), ]

## ploting
top_gm_genes <- dfg[order(-abs(dfg$avg_log2FC)), "Gene"][1:10]
top_gm_genes <- unique(c(top_gm_genes, txt_genes))

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
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
                        # boxedLabels = TRUE,
                        # selectLab = c("DEPP1", "BAG3"),   
                        selectLab = top_gm_genes,
                        colCustom = keyvals,
                        labSize = 5.5,
                        pointSize = 2,
                        
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
                        legendLabSize = 16,
                        
                        # thresholds
                        FCcutoff = 0.25,
                        pCutoff = 0.05,
                        ) +
    theme(legend.position=c(0.85, 0.2),
            axis.title=element_text(size=18),
            axis.line=element_line(size=1.5)) +
    guides(color = guide_legend(override.aes = list(size=3, alpha = 0.7)))  

png(paste0(vDir, "Fig1_GM_Volcano.png"), width = 600, height = 550)
print(p7)
dev.off()

