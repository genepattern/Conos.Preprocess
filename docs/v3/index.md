# Preprocess and Embed data with Conos

**Description**: Conos.Preprocess provides an entrypoint for the Clustering On Networks of Samples pipeline for integrating multiple single cell datasets into a common lower-dimensional embedding.

**Contact**: [genepattern.org/help](https://genepattern.org/help)

**Author(s)**:  Kharchenko Lab, Department of Biomedical Informatics, Harvard Medical School. <br> Wrapped as a module by Anthony S. Castanza, Edwin Juárez, and Barbara A. Hill Mesirov Lab, UCSD School of Medicine.

**Summary**: This is step 1 of 3 in the Conos pipeline. This module will read single cell data (e.g., count files) from different datasets and project each of those datasets into their own PCA space. The Conos.Preprocess module accepts either single cell gene expression matrices in tab delimited text with genes and rows and cells as columns, or Seurat datasets saved as RDS files. Tab delimited text matrices will be clustered and embedded into lower-dimensional space individually first using the Pagoda2 package, Seurat objects will have their embedding and clustering used as-is. The module then performs pairwise comparisons of the datasets to establish an initial error-prone mapping between cells of different datasets. These inter-sample edges are then combined with lower-weight intra-sample edges to construct a joint graph. This output can then be provided to the Conos.Cluster module to perform clustering on the global sample space.

**Source Publication**: Barkas N., Petukhov V., Nikolaeva D., Lozinsky Y., Demharter S., Khodosevich K., & Kharchenko P.V. 
Joint analysis of heterogeneous single-cell RNA-seq dataset collections. 
Nature Methods, (2019). doi:10.1038/s41592-019-0466-z <br>
See: [https://github.com/kharchenkolab/conos](https://github.com/kharchenkolab/conos) for full citation information

**Pipeline**: The complete pipeline utilizing the Conos workflow is available as a notebook at [Conos Integration for scRNA-seq](https://notebook.genepattern.org/hub/preview?id=442)

**Parameters:**

| Name | Description |
|:----------------|:----------------------------------|
| file list | Conos takes in text files (TXT), with tab-separated gene expression matrices (rows are genes, columns are unique cells), or preprocessed Seurat objects. Note that, while Conos can plot marker genes on embedded graphs from Seurat objects, it is not able to perform full differential expression analysis on Seurat data. Additionally, mixing of TXT datasets and Seurat objects in combined runs is not supported. |
| knn | Used by Pagoda2 if starting from raw matrices. Default number of neighbors to use in constructing the intra-sample kNN graph. |
| perplexity |  Used by Pagoda2 if starting from raw matrices. Perplexity to use in generating tSNE and largeVis embeddings (default=50) |
| pagoda odgenes | Used by Pagoda2 if starting from raw matrices. Number of top overdispersed genes to use for sample independent clustering (default=3e3)  |
| projection method | Projection method to use for performing pairwise sample comparisons supports PCA or CPCA |
| ncomps | How many (C)PCA components to use when performing the pairwise comparisons (default=50) |
| conos odgenes | Number of top overdispersed genes for Conos to use in performing pairwise comparisons and embeddings |


**Output**:

The output of the Conos.Preprocess module is an RDS object object that contains the main conos object under $con, the space that was used for the projection under $con_space (PCA or CPCA), and the data source under $data_source (seurat or matrix).

```
# Save an object to a file
saveRDS(list(con = con, con_space = con_space, data_source = data_source), "conos_preprocess_output.rds")
print("saved conos_preprocess_output.rds")
```

**Module Language**: R 

## Technical note:
This module uses the Conos Docker container vpetukhov/conos:version-1.4.4 wrapped in the container genepattern/conos:2.1 <br>
[The source code for this module is maintained in the Conos Dockerfile repository.](https://github.com/genepattern/docker-conos)
