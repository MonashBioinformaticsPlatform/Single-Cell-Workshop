# Single Cell RNA Seq Training Seminar

In partnership with 10X Genomics

Monday and Tuesday, 1st and 3rd of March

Instructors: Dr David Powell, Dr Paul Harrison, Dr Nick Wong, Dr Dee Deveson-Lucas, Nitika Kandhari, Peinan Zhao, Patrick Deng

## Day Two: Using Seurat (SingleR and Monocle)

This workshop, conducted by the Monash Bioinformatics Platform, will cover how to extend analysis to contemporary third-party tools, Seurat and Monocle. We will be walking through the 3K PBMC Dataset Tutorial and extend this to Monocle with some trajectory analysis.

Important links:

* [Installation and Setup instructions](installation.html)
* [PBMC Workshop material](pbmc3k_tutorial.html/)

* [Challenge solutions](solutions.html) (no peeking!)

### Recommended Computer Requirements:

System Requirements:

*Windows:*

* Windows 8.1 (64-bit) or later
* 4GB RAM
* SSD storage highly recommended
* Updated video/display drivers recommended

*macOS:*

* macOS 10.15 (Catalina) or later
* 4GB RAM
* SSD storage highly recommended

Install latest versions of:
* R
* RStudio
* Seurat
* SeuratWrappers
* Monocle

*Venue:* Virtual workshop conducted via Zoom with breakout rooms. Zoom links will be provided after registration is complete.



### Suggested Further Reading Material

* [Orchestrating Single Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/) - this book teaches single cell analysis with the bioconductor ecosystem of packages rather than Seurat. Regardless of your preference for Bioconductor or Seurat, it provides an excellent grounding and further depth and rationale behind each step of a single cell analysis.
* [Seurat tutorials for gene expression, spatial & multimodal analysis](https://satijalab.org/seurat/articles/get_started.html)
* [Getting started with Signac - the sibling package to Seurat for scATAC analysis](https://satijalab.org/signac/articles/overview.html)
* [Monocle documentationn for trajectories](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/)
* [Cell Annotation with SingleR](http://bioconductor.org/books/devel/SingleRBook/)
* [VDJ analysis with Immcantation](https://immcantation.readthedocs.io/en/stable/)

### Useful links arising from the discussion during the workshop
*	[10x Genomics link to ribosomal protein expression](https://kb.10xgenomics.com/hc/en-us/articles/218169723-What-fraction-of-reads-map-to-ribosomal-proteins-)
*	[10x Genomics link to mitochondrial gene expression](https://kb.10xgenomics.com/hc/en-us/articles/360001086611-Why-do-I-see-a-high-level-of-mitochondrial-gene-expression-)
*   [scRNA Tools, catelogue of tools for scRNA Seq analysis](https://www.scrna-tools.org/)
#### Data interpretation
*	[Interactive website explaining UMAP and comparision to t-SNE.](https://pair-code.github.io/understanding-umap/)
*	[OSCA, dimensionality reduction interpretation](http://bioconductor.org/books/3.14/OSCA.basic/dimensionality-reduction.html#visualization-interpretation)
#### Data tools and visualisation 
*	[scTransform Vignette](https://satijalab.org/seurat/articles/sctransform_vignette.html)
*	[Link to the workflowr library](https://github.com/jdblischak/workflowr)
*	[iSEE Bioconductor library, interactive explorer](https://bioconductor.org/packages/release/bioc/html/iSEE.html)
*	[ShinyCell makes interactive Shiny app from Seurat output](https://github.com/SGDDNB/ShinyCell)
*	[iCellR interactive data explorer](https://github.com/rezakj/iCellR)
*	[Diffusion maps for single cell instead of umaps](https://www.helmholtz-munich.de/icb/research/groups/marr-lab/software/destiny/index.html)
#### Papers
*   [Doublet cell detection method benchmarking paper.](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30459-2)
*	[From Louvain to Leiden: guaranteeing well-connected communities](https://www.nature.com/articles/s41598-019-41695-z)

#### Reference data and databases
*	[Gene tissue expression database](https://gtexportal.org/home/)
*	[ImmGen Database and Explorer](https://www.immgen.org/Databrowser19/DatabrowserPage.html)
*	[Single Cell Study Portal from The Broad](https://singlecell.broadinstitute.org/single_cell)
*	[Common ref data for cell indexing](http://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references)
*	[Azimuth is a Seurat-friendly reference-based annotation tool](https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC)
*	[Celaref, cell reference annotation tool](https://www.bioconductor.org/packages/release/bioc/html/celaref.html)


