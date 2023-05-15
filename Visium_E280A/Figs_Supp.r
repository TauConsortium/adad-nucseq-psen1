################################################################################
# Plots for Supp Figs

# conda activate py39
################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(cowplot)

# sub directories
cDir <- "/home/eger/projects/Brain_Visium/v3/clustering/"
vDir <- "/home/eger/projects/Brain_Visium/v3/plots/"
sDir <- "/home/eger/projects/Brain_Visium/v3/plots/Supp/"
dir.create(file.path(sDir), showWarnings = FALSE)
################################################################################
Neu <- c("RELN","VIP","GAD2","SYNPR","GAD1","CNR1","SYT1","SCG2","TAC3",
        "GABRG2","GABRA1","STMN2","DLX1","KCNC2","TMEM130","RAB3C","SST",
        "VSTM2A","SNAP25","ROBO2","CALB2","KIT","CNTNAP2","GABRB2","FSTL5",
        "NRXN3","SYT4","GRIA1","VSNL1","INA","NPY","GRIN2A","IGF1","PENK",
        "ELAVL2","MYT1L","KCNQ5","MEG3","NRIP3","CHGB","CLSTN2","SCN2A",
        "RAB3B","ZMAT4","NELL1","PNOC","ERBB4","SPHKAP","C11ORF87",
        "ADARB2","SLC4A10","KIAA1324","GRIN2B","BCL11A","CELF4","PNMA2",
        "DISP2","NYAP2","SV2B","SERPINI1","SLC2A13","RGS8","RTN1","NAP1L2",
        "CCK","C8ORF34","DYNC1I1","SRRM4","RBFOX1","SLC12A5","NDRG4",
        "ZNF804A","LPPR4","SLITRK4","GPR158","NDNF","KCNJ3","PCSK2",
        "CADPS","OLFM3","GABBR2","SULT4A1","GLRA2","SYT13","CACNA2D1",
        "GDA","SYNGR3","MAL2","PGM2L1","SLC7A14","GPR83","FRMPD4","NELL2",
        "RGS4","CSRNP3","DCX","ATP1A3","ST8SIA3","UCHL1","GAP43")

Oli <- c("PLP1","MOBP","CLDN11","MBP","UGT8","ERMN","MOG","MAG","OPALIN","CNP",
        "MAL","GPR37","TF","MYRF","GJB1","ASPA","ENPP2","BCAS1","LPAR1","FA2H",
        "ENPP6","APOD","CNTN2","CRYAB","KLK6","ERBB3","ANLN","SEPT4","PLEKHB1",
        "TMEFF2","ST18","PTGDS","PEX5L","SLAIN1","QDPR","PLLP","TMEM125","HHIP",
        "LGI3","TUBB4A","PLEKHH1","S1PR5","MAP6D1","GSN","EVI2A","EDIL3",
        "CMTM5","GJC3","CA14","NFASC","TPPP","TMEM88B","TRIM59","CDH19","APLP1",
        "NIPAL4","ADAMTS4","STMN4","S100B","CA2","PRR18","OLIG1","FOLH1",
        "NINJ2","NDRG1","SLC24A2","SGK2","GALNT6","KCNA1","SH3TC2","TTLL7",
        "SH3GL3","DOCK5","SCD","FEZ1","SLC44A1","RHOU","PPP1R16B","TSPAN2",
        "C10ORF90","TNFAIP6","NKAIN2","MOB3B","PRKCQ","PPP1R14A","PLA2G16",
        "DBNDD2","CDK18","PCDH9","ANO4","AGPAT4","OMG","FGFR2","TMEM63A",
        "GLTP","CCP110","PLEKHG3","RAB33A","PSAT1","ZNF536")
################################################################################

# List of 6 samples
samples <- c("C140_c12", "C226_c12", "C363_c12", "C364_c12",
             "C363_c1", "C364_c1")

# load objects
for (sample in samples) {
    load(paste0(cDir, sample, "/", sample, "_QC_cluster.Rdata"))
}

# integrated obj
load(paste0(cDir,"sample_integrated.Rdata"))

objs <- c(C140_c12, C363_c12, C364_c12, 
            C364_c1)

################################################################################
## Marker genes for 4 other samples ##

for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    # Spatial SCT - SNAP25 & MBP
    for (gene in c("SNAP25", "MBP")){
        p1 <- SpatialFeaturePlot(obj,
                                alpha =  c(0.3, 1), # min, max
                                image.alpha = 0.4,
                                crop = FALSE,
                                features = gene) + 
                                scale_fill_distiller(
                                    type = "seq", 
                                    palette = "Blues", 
                                    direction=1,
                                    guide = guide_colourbar(direction = "vertical",
                                                        frame.colour = "black",
                                                        frame.linewidth = 0.7)) +
                                theme(legend.position=c(0.09, 0.13), # x, y
                                    legend.text=element_text(size=15),
                                    legend.title=element_blank(),
                                    legend.background = element_rect(fill = NA))

        png(paste0(sDir, sample, "_", gene, ".png"))
        print(p1)
        dev.off()
    }
}

################################################################################
## Spatial plots for other 2 samples ##

for (obj in objs) {
    sample <- obj@meta.data$orig.ident[1]

    # change matter labels
    df <- obj@meta.data

    df$matter[df$matter == "GM"] <- "Gray Matter"
    df$matter[df$matter == "WM"] <- "White Matter"

    obj <- AddMetaData(obj, metadata = df)

    Idents(obj) <- as.factor(obj@meta.data$matter)

    p5 <- SpatialDimPlot(obj,
                    alpha = c(0.5, 1), # min, max
                    image.alpha = 0.4,
                    crop = FALSE,
                    cols = c("Gray Matter" = "#00C1AA", 
                                "White Matter" = "#FC61C9",
                                "Unknown" = "gray50")) +
                    theme(legend.position=c(0.2, 0.13), # x, y
                        legend.text=element_text(size=16),
                        legend.title=element_blank(),
                        legend.background = element_rect(fill = NA),
                        legend.key = element_blank()) +
                    scale_size_manual() +
                    guides(color = guide_legend(override.aes = list(
                        size = 10,
                        alpha = 1)))

    png(paste0(sDir, sample, "_spatial.png"))
    print(p5)
    dev.off()       
}

################################################################################
## UMAPs for clusters ##
sample_integrated@meta.data$seurat_clusters <- paste0("cluster ", as.character(sample_integrated@meta.data$seurat_clusters))
Idents(sample_integrated) <- as.factor(sample_integrated@meta.data$seurat_clusters)

p6 <- DimPlot(sample_integrated, 
            reduction = "umap",
            # cols = c("0" = "#00C1AA", # gray matter
            #         "1" = "#FC61C9", # white matter
            #         "2" = "#7997FF") # unknown
                    ) + 
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

png(paste0(sDir, "UMAP_integrated.png"))
print(p6)
dev.off()

################################################################################
## Integrated data - top cluster marker genes ##

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
de_markers <- FindAllMarkers(sample_integrated, 
                            only.pos = TRUE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)

# sort by logFC
df1 <- de_markers[order(-de_markers$avg_log2FC), ]

# remove cluster 2
df1 <- df1[!df1$cluster == 2, ]

# neuronal markers
top_Neu <- rownames(head(df1[rownames(df1) %in% Neu,], 10))

# oli markers
top_Oli <- rownames(head(df1[rownames(df1) %in% Oli,], 10))

# set assay
DefaultAssay(sample_integrated) <- "SCT"

# make feature plots
for (gene in top_Neu){
    p7 <- FeaturePlot(sample_integrated, 
                        features=gene) + 
                    NoLegend() +
                    theme(plot.title = element_text(
                            vjust = - 6,
                            size=22),
                        # legend.text=element_text(size=18),
                        axis.title.x = element_text(size=16),
                        axis.title.y = element_text(size=16),
                        axis.text.x=element_text(size=14),
                        axis.text.y=element_text(size=14),
                        axis.line = element_line(size = 1),
                        axis.ticks = element_line(size = 1)) +
                    scale_x_continuous(limits = c(-10, 10)) +
                    scale_y_continuous(limits = c(-5, 5),
                        breaks=c(-3, 0, 3))
    png(paste0(sDir, "UMAP_Neu_", gene, ".png"))
    print(p7)
    dev.off()
}

for (gene in top_Oli){
    p7 <- FeaturePlot(sample_integrated, 
                        features=gene) + 
                    NoLegend() +
                    theme(plot.title = element_text(
                            vjust = - 6,
                            size=22),
                        # legend.text=element_text(size=18),
                        axis.title.x = element_text(size=16),
                        axis.title.y = element_text(size=16),
                        axis.text.x=element_text(size=14),
                        axis.text.y=element_text(size=14),
                        axis.line = element_line(size = 1),
                        axis.ticks = element_line(size = 1)) +
                    scale_x_continuous(limits = c(-10, 10)) +
                    scale_y_continuous(limits = c(-5, 5),
                        breaks=c(-3, 0, 3))

    png(paste0(sDir, "UMAP_Oli_", gene, ".png"))
    print(p7)
    dev.off()
}


################################################################################
## Integrated data - Violin Plots ##

Idents(sample_integrated) <- as.factor(sample_integrated@meta.data$matter)

genes <- c("SNAP25", "CCK", "UCHL1", "SYT1")
for (gene in genes){
    p8 <- VlnPlot(sample_integrated, 
            idents = c("Gray Matter", "White Matter"),
            features=gene,
            cols = c("Gray Matter" = "#00C1AA", 
                    "White Matter" = "#FC61C9"),
            group.by = "matter",
            assay = "Spatial",
            slot = "data",
            pt.size = 0,
            y.max = 5) +
            labs(x = "", y = "log (counts)") + 
            NoLegend() +
            theme(plot.title = element_text(
                            vjust = -4,
                            size=22),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                axis.text.x=element_text(size=14),
                axis.text.y=element_text(size=14))
    
    png(paste0(sDir, "Violin_Neu_", gene, ".png"), width = 350, height = 500)
    print(p8)
    dev.off()
}

genes <- c("MBP", "MOBP", "TF", "CNP")
for (gene in genes){
    p8 <- VlnPlot(sample_integrated, 
            idents = c("Gray Matter", "White Matter"),
            features=gene,
            cols = c("Gray Matter" = "#00C1AA", 
                    "White Matter" = "#FC61C9"),
            group.by = "matter",
            assay = "Spatial",
            slot = "data",
            pt.size = 0,
            y.max = 7) +
            labs(x = "", y = "log (counts)") + 
            NoLegend() +
            theme(plot.title = element_text(
                            vjust = -4,
                            size=22),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                axis.text.x=element_text(size=14),
                axis.text.y=element_text(size=14))
    
    png(paste0(sDir, "Violin_Oli_", gene, ".png"), width = 350, height = 500)
    print(p8)
    dev.off()
}