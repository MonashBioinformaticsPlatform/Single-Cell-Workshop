
# This content is the [Seurat PBMC tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) with an additional sections on SingleR for cell type classification and Monocle for pseudotime trajectory analysis.

# Please go to the [installation page](https://monashbioinformaticsplatform.github.io/Single-Cell-Workshop/installation.html) for instructions on how to install the libraries used for this workshop. There are also instructions for downloading the [raw data](http://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) there as well.

# [The workshop homepage is here](https://monashbioinformaticsplatform.github.io/Single-Cell-Workshop/)


# Setup the Seurat Object --------

# For this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. The raw data can be found [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).

# We start by reading in the data. The Read10X() function reads in the output of the [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).

# We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For a technical discussion of the Seurat object structure, check out the [GitHub Wiki](https://github.com/satijalab/seurat/wiki). For example, the count matrix is stored in pbmc@assays$RNA@counts.

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#   **What does data in a count matrix look like?**

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D","TCL1A","MS4A1"), 1:30]

# The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0,  Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size / sparse.size


#### Discussion: The Seurat Object in R --------

# Lets take a look at the seurat object we have just created in R, pbmc

# To accomodate the complexity of data arising from a single cell RNA seq experiment, the seurat object keeps this as a container of multiple data tables that are linked.

# The functions in seurat can access parts of the data object for analysis and visualisation, we will cover this later on.

# There are a couple of concepts to discuss here.

# **Class**

# These are essentially data containers in R as a class, and can accessed as a variable in the R environment.

# Classes are pre-defined and can contain multiple data tables and metadata. For Seurat, there are three types.

# * Seurat - the main data class, contains all the data.
# * Assay - found within the Seurat object. Depending on the experiment a cell could have data on RNA, ATAC etc measured
# * DimReduc - for PCA and UMAP

# **Slots**

# Slots are parts within a class that contain specific data. These can be lists, data tables and vectors and can be accessed with conventional R methods.

# **Data Access**

# Many of the functions in Seurat operate on the data class and slots within them seamlessly. There maybe occasion to access these separately to hack them, however this is an advanced analysis method.

# The ways to access the slots can be through methods for the class (functions) or with standard R accessor nomenclature.

# **Examples of accessing a Seurat object**

# The assays slot in pbmc can be accessed with pbmc@assays.

# The RNA assay can be accessed from this with pbmc@assays$RNA.

# We often want to access assays, so Seurat nicely gives us a shortcut pbmc$RNA. You may sometimes see an alternative notation pbmc[["RNA"]].

# In general, slots that are always in an object are accessed with @ and things that may be different in different data sets are accessed with $.

# **Have a go**

# Use str to look at the structure of the Seurat object pbmc.

# What is in the meta.data slot within your Seurat object currently? What type of data is contained here?

# Where is our count data within the Seurat object?


# Standard pre-processing workflow --------

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.


## QC and selecting cells for further analysis --------

# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics [commonly used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community include

# * The number of unique genes detected in each cell.
#     + Low-quality cells or empty droplets will often have very few genes
#     + Cell doublets or multiplets may exhibit an aberrantly high gene count
# * Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# * The percentage of reads that map to the mitochondrial genome
#     + Low-quality / dying cells often exhibit extensive mitochondrial contamination
#     + We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
#     + We use the set of all genes starting with MT- as a set of mitochondrial genes

# The $ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")


#### Challenge: The meta.data slot in the Seurat object --------

# Where are QC metrics stored in Seurat?

# * The number of unique genes and total molecules are automatically calculated during CreateSeuratObject()
#     + You can find them stored in the object meta data

# What do you notice has changed within the meta.data table now that we have calculated mitochondrial gene proportion?

# Could we add more data into the meta.data table?

# ###

# In the example below, we visualize QC metrics, and use these to filter cells.

# * We filter cells that have unique feature counts over 2,500 or less than 200
# * We filter cells that have >5% mitochondrial counts


#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Lets look at the number of features (genes) to the percent mitochondrial genes plot.

plot3 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot3


#### Challenge: Ribosomal gene expression as a QC metric --------

# Ribosomal gene expression could be another factor to look into your cells within your experiment.

# Create more columns of metadata using PercentageFeatureSet function, this time search for ribosomal genes. We can  calculate the percentage for the large subunit (RPL) and small subunit (RPS) ribosomal genes.

# Use FeatureScatter to plot combinations of metrics available in metadata. How is the mitochondrial gene percentage related to the ribosomal gene percentage? What can you see? Discuss in break out.

# **Code for challenge**
# Create new meta.data columns to contain percentages of the large and small ribosomal genes.

# Then plot a scatter plot with this new data. You should find that the large and small ribosomal subunit genes are correlated within cell.

# What about with mitochondria and gene, feature counts?

# These are the cells you may want to exclude.

# **Advanced Challenge**
# Highlight cells with very low percentage of ribosomal genes, create a new column in the meta.data table and with FeatureScatter make a plot of the RNA count and mitochondrial percentage with the cells with very low ribosomal gene perentage.

# ###

# Okay we are happy with our thresholds for mitochondrial percentage in cells, lets apply them and subset our data. This will remove the cells we think are of poor quality.

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Lets replot the feature scatters and see what they look like.

plot5 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot5 + plot6


# Normalizing the data --------

# After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc$RNA@data.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)

# For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. However, this isn't required and the same behavior can be achieved with:

pbmc <- NormalizeData(pbmc)


# Identification of highly variable features (feature selection) --------

# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and [others](https://www.nature.com/articles/nmeth.2645) have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

# Our procedure in Seurat is described in detail [here](https://doi.org/10.1016/j.cell.2019.05.031), and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#### Challenge: Labelling Genes of Interest --------

# What if we wanted to look at genes we are specifically interested in? We can create a character vector of gene names and apply that to this plot.

# Make a plot with labels for the genes IL8, IDH2 and CXCL3.


# Scaling the data --------

# Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:

# * Shifts the expression of each gene, so that the mean expression across cells is 0
# * Scales the expression of each gene, so that the variance across cells is 1
#     + This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# * The results of this are stored in pbmc$RNA@scale.data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#   **This step takes too long! Can I make it faster?**

# Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, the default in ScaleData() is only to perform scaling on the previously identified variable features (2,000 by default). To do this, omit the features argument in the previous function call, i.e.

# pbmc <- ScaleData(pbmc)

# Your PCA and clustering results will be unaffected. However, Seurat heatmaps (produced as shown below with DoHeatmap()) require genes in the heatmap to be scaled, to make sure highly-expressed genes don't dominate the heatmap. To make sure we don't leave any genes out of the heatmap later, we are scaling all genes in this tutorial.

#   **How can I remove unwanted sources of variation, as in Seurat v2?**

# In Seurat v2 we also use the ScaleData() function to remove unwanted sources of variation from a single-cell dataset. For example, we could 'regress out' heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination. These features are still supported in ScaleData() in Seurat v3, i.e.:

# pbmc <- ScaleData(pbmc, vars.to.regress = 'percent.mt')

# However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow, SCTransform(). The method is described in our [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1), with a separate vignette using Seurat v3 [here](sctransform_vignette.html). As with ScaleData(), the function SCTransform() also includes a vars.to.regress parameter.


# Perform linear dimensional reduction --------

# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(pbmc$pca, dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')

# In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the 'dimensionality' of the dataset --------

# To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset. However, how many components should we choose to include? 10? 20? 100?

# In [Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8), we implemented a resampling test inspired by the JackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores, and repeat this procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value features.

# Note:
# The Seurat defaults are num.replicates=100, prop.freq=0.01.
# The parameters we're using here are just for speed,
#   and will give conservative p-values.

pbmc <- JackStraw(pbmc, num.replicate=10, prop.freq=0.1)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 'Significant' PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.

JackStrawPlot(pbmc, dims = 1:15)

# Seurat plots this in a weird way. Let's adjust it a bit.
JackStrawPlot(pbmc, dims = 1:15) +
  coord_cartesian() +
  geom_abline(intercept=0, slope=0.05)

# For each PC, the furthest point to the right below the solid line
#   gives proportion of significant genes with FDR 0.05.

# An alternative heuristic method generates an 'Elbow plot': a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an 'elbow' around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

ElbowPlot(pbmc)

# Identifying the true dimensionality of a dataset -- can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-12 as a cutoff.

# We chose 10 here, but encourage users to consider the following:

# * Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
# * We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). As you will observe, the results often do not differ dramatically.
# * We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results.


# Cluster the cells --------

# Seurat v3 applies a graph-based clustering approach, building upon initial strategies in ([Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8)). Importantly, the *distance metric* which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partitioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [[SNN-Cliq, Xu and Su, Bioinformatics, 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract) and CyTOF data [[PhenoGraph, Levine *et al*., Cell, 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'.

# As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors() function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

# To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [[SLM, Blondel *et al*., Journal of Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters() function implements this procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents() function.

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# Run non-linear dimensional reduction (UMAP/tSNE) --------

# Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = "umap-learn")
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = 'umap')

# You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.

saveRDS(pbmc, file = "pbmc_tutorial.rds")


#### Challenge: Try different cluster settings --------

# Run FindNeighbours and FindClusters again, with a different number of dimensions or with a different resolution. Examine the resulting clusters using DimPlot.

# To maintain the flow of this tutorial, please put the output of this exploration in a different variable, such as pbmc2!


# Finding differentially expressed features (cluster biomarkers) --------

# Seurat can help you find markers that define clusters via differential expression. By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells.  FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

# Seurat has several tests for differential expression which can be set with the test.use parameter (see our [DE vignette](de_vignette.html) for details). For example, the ROC test returns the 'classification power' abs(AUC-0.5)*2 for any individual marker, ranging from 0 = random to 1 = perfect.

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression. VlnPlot() (shows expression probability distributions across clusters), and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

#   **Other useful plots**
# These are ridgeplots, cell scatter plots and dotplots. Replace FeaturePlot with the other functions.

RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# For CellScatter plots, will need the cell id of the cells you want to look at. You can get this from the cell metadata (pbmc@meta.data).

head( pbmc@meta.data )
CellScatter(pbmc, cell1 = "AAACATACAACCAC-1", cell2 = "AAACATTGAGCTAC-1")

# DotPlots

DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Which plots do you prefer? Discuss.

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


# Assigning cell type identity to clusters --------

# Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:

# Cluster ID | Markers       | Cell Type
# -----------|---------------|----------
# 0          | IL7R, CCR7    | Naive CD4+ T
# 1          | CD14, LYZ     | CD14+ Mono
# 2          | IL7R, S100A4  | Memory CD4+
# 3          | MS4A1         | B
# 4          | CD8A          | CD8+ T
# 5          | FCGR3A, MS4A7 | FCGR3A+ Mono
# 6          | GNLY, NKG7    | NK
# 7          | FCER1A, CST3  | DC
# 8          | PPBP          | Platelet

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "pbmc3k_final.rds")

write.csv(x = t(as.data.frame(all_times)), file = "pbmc3k_tutorial_times.csv")


# SingleR --------

#install.packages("BiocManager")
#BiocManager::install(c("SingleCellExperiment","SingleR","celldex"),ask=F)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

# In this workshop we have focused on the Seurat package.  However, there is another whole ecosystem of R packages for single cell analysis within Bioconductor.  We won't go into any detail on these packages in this workshop, but there is good material describing the object type online : [OSCA](https://robertamezquita.github.io/orchestratingSingleCellAnalysis/data-infrastructure.html).

# For now, we'll just convert our Seurat object into an object called SingleCellExperiment.  Some popular packages from Bioconductor that work with this type are Slingshot, Scran, Scater.

sce <- as.SingleCellExperiment(pbmc)

# We will now use a package called SingleR to label each cell.  SingleR uses a reference data set of cell types with expression data to infer the best label for each cell.  A convenient collection of cell type reference is in the celldex package which currently contains the follow sets:

ls('package:celldex')

# In this example, we'll use the HumanPrimaryCellAtlasData set, which contains high-level, and fine-grained label types. Lets download the reference dataset

ref.set <- celldex::HumanPrimaryCellAtlasData()
head(unique(ref.set$label.main))

# An example of the types of "fine" labels.

head(unique(ref.set$label.fine))

# Now we'll label our cells using the SingleCellExperiment object, with the above reference set.

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)

# Keep any types that have more than 10 cells to the label, and put those labels back on our Seurat object and plot our on our umap.

lbls.keep <- table(pred.cnts$first.labels)>10
pbmc$SingleR.labels <- ifelse(lbls.keep[pred.cnts$first.labels], pred.cnts$first.labels, 'Other')
DimPlot(pbmc, reduction='umap', group.by='SingleR.labels')

# It is nice to see that SingleR does not use the clusters we computed earlier, but the labels do seem to match those clusters reasonably well.


# Using Monocle For Pseudotime Trajectory (Time permits) --------

# For this workshop, we'll use the PBMC data object with Monocle for pseudotime trajectory analysis. It's debatable whether this is a suitable dataset but will suit our needs for demonstration purposes.

# This content is based off the [Calculating Trajectories with Monocle 3 and Seurat](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html) material as well the [Monocle3 documentation](https://cole-trapnell-lab.github.io/monocle3/docs/starting/), combining it with the PBMC dataset from the original Seurat vignette. We recommend reading the [Monocle3 documentation](https://cole-trapnell-lab.github.io/monocle3/docs/starting/) for greater understanding of the Monocle package.

# Firstly, load Monocle:

library(monocle3)
library(SeuratWrappers)

# As the PBMC data has been processed, we can proceed with converting the pbmc Seurat object to a cell_data_set object, which is a class from the Monocle package. The as.cell_data_set function is used from the SeuratWrappers library and is used to convert the Seurat object into a cell_data_set object.

# While we have performed the general analysis steps of quality control, scaling and normalization, dimensionality reduction and clustering with Seurat, Monocle is also capable of performing these steps with its own in-built functions. It is often a matter of preference which package to use, depending on what downstream tasks the analyst would like to perform.

# ![](https://cole-trapnell-lab.github.io/monocle3/images/monocle3_new_workflow.png)
# We aren't going to delve deeply into the properties of the cell_data_set object. Just be aware that this is a different way to represent the count assay data and dimensionality reduction data. The functions from the Monocle package expects the scRNA data to be this class and therefore, the Seurat object needs to be converted to this class. It also means that the Seurat functions that we've been using will not work with the cell_data_set object.

cds <- as.cell_data_set(pbmc)

# While we have previously clustered the pbmc dataset using Seurat, Monocle will also calculate 'partitions' - these are superclusters of the Louvain/Leiden communties that are found using a kNN pruning method. The warning message during the conversion notes that Seurat doesn't calculate partitions and clusters need to be re-calculated using Monocle.

# Examine the cds object:

## Inspect the cds object and compare it to the Seurat pbmc object
cds

# Now re-cluster the cds object:

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# We can see that in the first plot, Monocle has identified 3 clusters but all the clusters fall within the same partition. Ideally, partitions should correspond to clusters of cells within the same path of differentiation or cell within the same trajectory.

# ![](https://ars.els-cdn.com/content/image/1-s2.0-S0263931913000495-gr1.jpg)

# *Source: [Haematopoiesis and red blood cells](https://www.sciencedirect.com/science/article/pii/S0263931913000495)*

# We can see from this figure of haematopoiesis that our PBMC sample contains a mix of cells from different cell types and are unlikely to be suitable for calculating a pseudotime trajectory. Nonetheless, we'll demonstrate the steps involved.

# Next, we need to run learn_graph to learn the trajectory graph. This function aims to learn how cells transition through a biological program of gene expression changes in an experiment.

cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

# As expected from the partition plot, Monocle thinks all the cells are from the same partition and therefore has plotted a trajectory line that connects all clusters.

# Can we fix this?

# Monocle currently thinks that all cells belong to the same partition. We might be able to tweak the clustering for a better result. One thing we can think about is that we have a different number of clusters generated by Monocle (3) when Seurat gave us 9.

# If we examine the default parameters by ?FindClusters and ?cluster_cells, we might notice that Seurat's default clustering algorithm is louvain while Monocle's is leiden. We aren't going to delve into the details of these algorithms, but we will say, just be aware of the default behavior of your analysis tools and that the choice in algorithm will affect the results of the clustering.

# We can change algorithm with cds <- cluster_cells(cds, cluster_method = "louvain") but in this case, we might just try altering the resolution with the default leiden algorithm to increase the number of clusters yielded. Changing the k argument will change the number of nearest neighbors used when creating the k nearest neighbor graph. A large k value (the default is 20) reduces the number of clusters (therefore the bigger k is, the less clusters will be generated) and vice versa (smaller k value - more clusters).

cds <- cluster_cells(cds, k = 5, random_seed = 5)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# The UMAP on the left looks under clustered compared to our original clustering with Seurat. We'd probably need to tweak more parameters to get Monocle to match the Seurat clustering. We don't necessarily need to do that because our cds object actually still has the meta-data about the Seurat clusters stored in it (examine this with head(colData(cds)). However, importantly, our partition plot looks a little more sensible and no longer has lumped all cells into one supercluster.

# Let's re-run the learn_graph step:

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)

# Monocle now 'correctly' builds trajectories that recognizes distinct cell lineages.

# We might choose to remove the B-cells and monocytes and focus just on the cluster of CD4T/CD8T cells, as this is the largest group of cells.

# Create a vector of idents to keep
selected_ids <- c("Naive CD4 T", "Memory CD4 T", "CD8 T")
tcells_pbmc <- subset(pbmc, idents = selected_ids ) ## subset the PBMC seurat object to tcells
cds <- as.cell_data_set(tcells_pbmc) ## convert this to cell_data_set
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)

# The next step is to order cells in pseudotime:

# >Pseudotime is a measure of how much progress an individual cell has made through a process such as cell differentiation.
# >
# > In many biological processes, cells do not progress in perfect synchrony. In single-cell expression studies of processes such as cell differentiation, captured cells might be widely distributed in terms of progress. That is, in a population of cells captured at exactly the same time, some cells might be far along, while others might not yet even have begun the process. This asynchrony creates major problems when you want to understand the sequence of regulatory changes that occur as cells transition from one state to the next. Tracking the expression across cells captured at the same time produces a very compressed sense of a gene's kinetics, and the apparent variability of that gene's expression will be very high.
# >
# > By ordering each cell according to its progress along a learned trajectory, Monocle alleviates the problems that arise due to asynchrony. Instead of tracking changes in expression as a function of time, Monocle tracks changes as a function of progress along the trajectory, which we term "pseudotime". Pseudotime is an abstract unit of progress: it's simply the distance between a cell and the start of the trajectory, measured along the shortest path. The trajectory's total length is defined in terms of the total amount of transcriptional change that a cell undergoes as it moves from the starting state to the end state.

# *Source: [Monocle's documentation](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/#order-cells)*

# Monocle needs to be told where the 'beginning' of the biological process is. There are a variety of ways that this can be determined - the Monocle documentation has a custom function to find the root of the trajectory based on a subset of cells. If the order_cells function is used without providing which cells to use, it will launch an interface in which we can directly select cells we think are at the beginning of the trajectory.

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds,  cell_type="Naive CD4 T"){
  cell_ids <- which(colData(cds)[, "ident"] == cell_type)

  closest_vertex <-
  cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)$UMAP)$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# We can now plot the trajectory and color cells by pseudotime:

plot_cells(cds, color_cells_by = "pseudotime",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "ident",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE)

# Discuss the results of this pseudotime trajectory (remembering this is a bogus example):

# * CD8T cells are not the furthest cell type from the naive CD4 t-cells
# * Naive CD4 T-cells get split into two groups - cells at the root state and 'early' in terms of pseudotime and then cells that are at the end of the pseudotime  timeline
# * Would you interpret this as naive CD4 T cells shifting into memory CD4 T-cells then CD8T cells and then back to CD4 naive?
# * An analysis tool will always try to give you some sort of answer - it's important to think about whether the tool we're using is appropriate for a given dataset
# * What happens if we left in the NK cells?

#   **Session Info**

sessionInfo()
