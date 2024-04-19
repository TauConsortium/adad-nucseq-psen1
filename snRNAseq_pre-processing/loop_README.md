# Pre-Processing Pipeline for snRNA Sequencing
#### Developed by Davis Westover based on code by Camila Almeidas, meant to be run in R

### Summary 
 This code handles pre-processing for snRNA data. Samples are processed separately before they are able to be merged into one dataset. This involves subsetting out low quality cells, normalizing data, and subsetting doublets.

### Software Needed

 - Download R and Rstudio
 - Download Cellranger or use the 10x genomics cloud
 - RPackages dplyr, Seurat, Matrix, ggplot2, sctransform, EnhancedVolcano, DoubletFinder, and pheatmap
 - Doublet finder can be downloaded with remotes from https://github.com/chris-mcginnis-ucsf/DoubletFinder

### Useful guides / Acknowledgements
-[Seurat v4 tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial/)\
-[doubletfinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder/)\
-[sctransform](https://satijalab.org/seurat/articles/sctransform_vignette.html/)

### Process Overview:
A) **Subset out low quality cells** \
B) **Normalize data** \
C) **Doublet finder** \
D) **Subset out doublets** \
E) **Merge samples** 


### Data needed:
- You will need already sequenced snRNA Data in a fastq file. 
  - Single nucleus RNA sequencing is used to profile gene expression. Cell walls of nuclei are lysed (ruptured) to capture RNA, which can be sorted and sequenced too help profile the genes expressed in the cells we took nuclei from. This sequencing was done for us by a lab on samples from 8 carriers of the PSEN1-E280A gene, as well as 8 non mutation carriers with sporadic AD and controls with a range of APOE genotypes.
- Cellranger will generate feature-barcode matrices, determine clusters, and perform gene expression analysis.
    - This compiles the RNA data collected and marks cells with the information collected, as well as clusters what the cells likely are.
- Create a directory to store the resulting file.

---
### Variables Used

| **Variables** |                                                    **Description**                                                   |
|:----------------:|:--------------------------------------------------------------------------------------------------------------------:|
|       percent.mt    | Proportion of transcripts mapping to mitochondrial genes |
|        nFeature_RNA  | Count of unique genes in each cell |
|        nCount_RNA      | Count of genes in each cell  |
|        homotypic.prop    | Leverages cell annotations to model the proportion of homotypic doublets, modeled as the sum of squared annotation frequencies  |
|        nExp_poi    | Defines the pANN threshold used to make final doublet/singlet predictions |
|        nExp_poi.adj  | Adjusted nExp_poi, calculated with nExp_poi and homotypic.prop |
|        pANN_col_name  | Calculated pANN (proportional of artificial k nearest neighbors) value |
|        MT.genes    | List of mitochondrial genes  |
___


### Steps in R

**1. Create the Seurat objects**
- Using the Seurat package, read in your files with Read10x.
- Convert these objects into Seurat Objects.
- Initial pre-processing closely follows Seurat's clustering guide

**2. Visualize and subset low quality cells**
  - 'Low quality cells' often have very few genes, and few unique genes. These genes often exhibit extensive mitochondrial contamination. We use these parameters to remove a subset of the cells.
  - Doublets can also be sought out by noticing genes with an aberrantly high gene count. Doublets are cells counted twice and must be filtered out the dataset.
  - The graphs should roughly follow this improvement before and after filtering (note the x-axis changes):

![pre_and_post_feature_counts](https://github.com/daviswestover/read_better/assets/91497472/ecc50afb-d13f-454c-8443-14e654f02f27)

**3. Run SCTransform, PCA, and UMAP** 
- SCTransform does three things: 
    - Normalize Data: Divides each count by total cells, multiplies by a scale factor, then log-transforms data. This adjusts the scale of the data to make it easier to work with, while still proportionally accounting for the numeric differences. 
    - Scaling Data: Centers and scales the data, further improving numeric stability. It is important to do this before PCA.
    - FindVariableFeatures: Identifies variable features that have strong affects on variation, while controling for mean variability and average expression.
- RunPCA() (Primary Component Analysis) determines the primary genes having an effect on cell results.
- RunUmap() (Uniform Manifold Approximation and Projection) is a dimnesion reduction technique.
    
**4. Run PK identification**
  - PK defines Primary Component neighborhood size used to compute pANN (Proportion of artificial nearest neighbors.) 
  - PK value for our samples followed the doublet finder guide's pk value instead of individually identifying values. However, bcmvn_sample can be used to select a custom pk value.
  
**5. Homotypic Doublet Proportion Estimate**
  - nExp defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities already collected from cellbender.

**6. Run doublet finder**
- It is important that low quality cells are removed before running doublet finder.
- The first run establishes a pann value which is used when run the second time.
    
**7. Subset your Seurat Object to exclude doublet cells and MTgenes**
- A doubletfinder column is created after running doubletfinder, we mark this as df_class_col_name and use it to identify marked doublets.
- Identified cells are subsetted out.
- Mitochondrial genes are not relevant to the work and can be taken out.
    - Create a list of all genes that include 'MT', along with the gene S402, which is also unnecessary.
    - Subset sample genes to disclude these genes.

**8. Merge the separately processed data into one Seurat object**
- With the data now all processed separately, use merge to create one object of all the samples.
- Save this work with saveRDS. Data is now ready for analysis.

