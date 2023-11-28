# Single nucleus RNA sequencing & Spatial Transcriptomics analyses from brain tissue
Repository containing the computational code for the "Single Nucleus RNA Sequencing Demonstrates an Autosomal Dominant Alzheimer’s Disease Profile and Possible Mechanism of Disease Protection" manuscript by Almeida et al.

Computational pipeline developed by:\
Camila Almeida <caalmeida@ucsb.edu>\
Sarah Eger <eger@ucsb.edu>\
Caroline He <caroline_he@berkeley.edu>\
Davis Westover <daviswestover@ucsb.edu>

**Quality Control and Cell clustering**

[snRNAseq.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/snRNAseq.R) \

**Differential Gene Expression analyses**

[APOE3_Christchurch_GLMM/1_glmer.nb_all_cts_prep.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/APOE3_Christchurch_GLMM/1_glmer.nb_all_cts_prep.r)
Takes a seurat object as input and outputs a dataframe per cell type that includes the necessary metadata and gene counts for every cell to be included in the negative-binomial generalized linear mixed model. This was run twice (1) with 8 individuals (1 Christchurch homozygote & 7 E280A Christchurch non-carriers) and (2) with 10 individuals (3 Christchurch heterozygotes & 7 E280A Christchurch non-carriers). Metadata variable selection was adjusted accordingly (see section ## Formatting metadata).

[APOE3_Christchurch_GLMM/2_glmer.nb_single_ct_parallelized.r](https://github.com/acostauribe/nucseq-PSEN1-E280A/blob/main/APOE3_Christchurch_GLMM/2_glmer.nb_single_ct_parallelized.r)
Takes a single dataframe for a given cell type and outputs the results of the negative-binomial generalized linear mixed model run on every gene. Multiple genes are run in parallel to speed up the process. Change the number of cores depending on your machine. This was run separately for every cell type for each group comparison.


**High dimensional weighted gene co-expression network analysis (hd-WGCNA)**

[hdWGCNA/hdWGCNA-module_identification.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-module_identification.R)
This code is run after basic seurat object quality control and preprocessing. First, groups of individual cells within the Seurat object are aggregated into 'metacells', and this aggregated data is used to identify correlations between genes across metacells of the same celltype. Groups of highly interconnected genes are called 'modules'. After modules are identified, the top hub genes for each module is visualized, and statistical tests are run to identify modules that are differentially expressed between different conditions.

[hdWGCNA/hdWGCNA-Ast_marker_overlap.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-Ast_marker_overlap.R)
This code is run after hdWGCNA-module_identification.R. At this point, each celltype has particular modules associated with it. Each celltype also is clustered into subclusters. In the case of the astrocytes, they have 10 modules associated with it and the astrocyte population can be subsetted into 4 different clusters. In this script, the module gene lists, and the marker lists for each of the clusters, are compared. The significance of the list overlap is visualized, and we can identify particularly significant overlaps.

**Spatial Transcriptomics**
