# Pre-Processing Pipeline for snRNA Sequencing
#### Developed by Davis Westover based on code by Camila Almeida

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
A **Subset out low quality cells**
B **Normalize data**
C **Doublet finder**
D **Subset out doublets**
E **Merge samples**


### Data needed:
- You will need already sequenced snRNA Data in a fastq file.
- Cellranger will generate feature-barcode matrices, determine clusters, and perform gene expression analysis.
- Create a directory to store the resulting file

### Steps in R
Now, using the Seurat package, Read in your file with Read10x.
Initial pre-processing closely follows Seurat's clustering guide:

**1. Create the Seurat objects**

**2. Visualize and subset low quality cells**
  -'Low quality cells' often have very few genes, and few unique genes. These genes often exhibit extensive mitochondrial contamination
  - Doublets can also be sought out be noticing genes with an aberrantly high gene count. The expected number of doubles is 1 for every thousand cells total.
  
  
**3. Run SCTransform, PCA, and UMAP**
  - 'Variable features' subsets features that exhibit high cell-to-cell variation in the dataset
    
**4. Run pk identification**
  - PK value for our samples followed the doublet finder guide's pk value instead of individually identifying values

**5. Run doublet finder**
    -Important low quality cells are removed before doublet finder
    -First time running establishes pann value which is used when run the second time
    
**6. Subset your Seurat Object to exclude doublet cells**

**7. Subset to remove MTgenes (and gene S402) from your Seurat object**
