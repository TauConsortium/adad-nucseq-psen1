# Single cell nuclear RNA analyses from brain tissue
Repository containing the computational code for the "Chaperone Gene Expression Distinguishes E280A PSEN1 Autosomal Dominant Alzheimer´s by Single Nucleus RNA Sequencing" manuscript by Almeida et al.

Computational pipeline developed by:\
Camila Almeida <caalmeida@ucsb.edu>\
Sarah Eger <eger@ucsb.edu>\
Caroline He <caroline_he@berkeley.edu>\
Juliana Acosta-Uribe <acostauribe@ucsb.edu>

**Quality Control and Cell clustering**

[snRNAseq.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/snRNAseq.R) \
Camila please descirbe your code here

**Differential Gene Expression analyses**

Pseudo bulk
linear mixed model (2.2 code)

**High dimensional weighted gene co-expression network analysis (hd-WGCNA)**

[hdWGCNA-module_identification.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-module_identification.R)
This code is run after basic seurat object quality control and preprocessing. First, groups of individual cells within the Seurat object are aggregated into 'metacells', and this aggregated data is used to identify correlations between genes across metacells of the same celltype. Groups of highly interconnected genes are called 'modules'. After modules are identified, the top hub genes for each module is visualized, and statistical tests are run to identify modules that are differentially expressed between different conditions.

[hdWGCNA-Ast_marker_overlap.R](https://github.com/KosikLabUCSB/nucseq-PSEN1-E280A/blob/main/hdWGCNA/hdWGCNA-Ast_marker_overlap.R)
This code is run after hdWGCNA-module_identification.R. At this point, each celltype has particular modules associated with it. Each celltype also is clustered into subclusters. In the case of the astrocytes, they have 10 modules associated with it and the astrocyte population can be subsetted into 4 different clusters. In this script, the module gene lists, and the marker lists for each of the clusters, are compared. The significance of the list overlap is visualized, and we can identify particularly significant overlaps.

**Spatial Transcriptomics**
