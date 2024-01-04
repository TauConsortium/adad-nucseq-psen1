# Single nucleus & Spatial transcriptomics analyses

Repository containing the computational code for the "Single Nucleus RNA Sequencing Demonstrates an Autosomal Dominant Alzheimer’s Disease Profile and Possible Mechanism of Disease Protection" manuscript by Almeida et al.

https://api.github.com/repos/TauConsortium/adad-nucseq-psen1
[![DOI](https://zenodo.org/badge/{github_id}.svg)](https://zenodo.org/badge/latestdoi/{github_id})

**Computational pipeline developed by:**\
Camila Almeida <caalmeida@ucsb.edu>\
Sarah J. Eger <eger@ucsb.edu>\
Caroline He <caroline_he@berkeley.edu>\
Davis Westover <daviswestover@ucsb.edu>\
Juliana Acosta-Uribe <acostauribe@ucsb.edu>

### **Quality Control for expression matrix**

[snRNAseq_pre-processing/pre-processing_loop.R](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/snRNAseq_pre-processing/pre-processing_loop.R) has separate readme.md: [snRNAseq_pre-processing/loop_README.md](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/snRNAseq_pre-processing/loop_README.md)

### **Data Processing, Analyses, Visualization, and Differential Expression Testing**

[snRNAseq.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/snRNAseq.R)
Data normalization, scaling, integration, clustering, cell-type identification.

### **Identification of differentially expressed genes in cell-type subpopulations**

[APOE3_Christchurch_GLMM/1_glmer.nb_all_cts_prep.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/APOE3_Christchurch_GLMM/1_glmer.nb_all_cts_prep.r)
Takes a seurat object as input and outputs a dataframe per cell type that includes the necessary metadata and gene counts for every cell to be included in the negative-binomial generalized linear mixed model. This was run twice (1) with 8 individuals (1 Christchurch homozygote & 7 E280A Christchurch non-carriers) and (2) with 10 individuals (3 Christchurch heterozygotes & 7 E280A Christchurch non-carriers). Metadata variable selection was adjusted accordingly (see section ## Formatting metadata).

[APOE3_Christchurch_GLMM/2_glmer.nb_single_ct_parallelized.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/APOE3_Christchurch_GLMM/2_glmer.nb_single_ct_parallelized.r)
Takes a single dataframe for a given cell type and outputs the results of the negative-binomial generalized linear mixed model run on every gene. Multiple genes are run in parallel to speed up the process. Change the number of cores depending on your machine. This was run separately for every cell type for each group comparison.

### **High dimensional weighted gene co-expression network analysis (hd-WGCNA)**

[hdWGCNA/hdWGCNA-module_identification.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-module_identification.R)
This code is run after basic seurat object quality control and preprocessing. First, groups of individual cells within the Seurat object are aggregated into 'metacells', and this aggregated data is used to identify correlations between genes across metacells of the same celltype. Groups of highly interconnected genes are called 'modules'. After modules are identified, the top hub genes for each module is visualized, and statistical tests are run to identify modules that are differentially expressed between different conditions.

[hdWGCNA/hdWGCNA-Ast_marker_overlap.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-Ast_marker_overlap.R)
This code is run after hdWGCNA-module_identification.R. At this point, each celltype has particular modules associated with it. Each celltype also is clustered into subclusters. In the case of the astrocytes, they have 10 modules associated with it and the astrocyte population can be subsetted into 4 different clusters. In this script, the module gene lists, and the marker lists for each of the clusters, are compared. The significance of the list overlap is visualized, and we can identify particularly significant overlaps.

### **Spatial Transcriptomics**

[Spatial_Transcriptomics/1_Mapping_spatial_gene_expression_reads_SPACERANGER.py](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/1_Mapping_spatial_gene_expression_reads_SPACERANGER.py)
Mapping spatial gene expression reads to the genome and microscope images with spaceranger.

[Spatial_Transcriptomics/2_compare_spaceranger_auto_and_manual_alignment.py](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/2_compare_spaceranger_auto_and_manual_alignment.py) Compare auto and manual alignment using the metrics outputed by spaceranger.

[Spatial_Transcriptomics/3_preprocess_individual_samples_SEURAT.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/3_preprocess_individual_samples_SEURAT.r) Preprocess and cluster each sample individually with Seurat.

[Spatial_Transcriptomics/4_integrate_and_cluster_SEURAT.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/4_integrate_and_cluster_SEURAT.r) Integrate spots from all samples and cluster with Seurat.

[Spatial_Transcriptomics/5a_glmer.nb_WM_GM_prep.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/5a_glmer.nb_WM_GM_prep.r) (a) Prepare dataframes
[Spatial_Transcriptomics/5b_glmer.nb_WM_GM_parallelized_LME4.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/5b_glmer.nb_WM_GM_parallelized_LME4.r) (b) run a negative binomial GLMM on each gene
[Spatial_Transcriptomics/5c_glmer.nb_4_sample_LME4.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/5c_glmer.nb_4_sample_LME4.r) (c) run without outlier sample.

[Spatial_Transcriptomics/6_filter_significant_DEGs.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/6_filter_significant_DEGs.r) Filter significant differentially expressed genes.

[Spatial_Transcriptomics/7_fisher_test_DEGs.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/Spatial_Transcriptomics/7_fisher_test_DEGs.r) Perform Fisher's Exact Test to check for significant overlap between lists of differentially expressed genes.
