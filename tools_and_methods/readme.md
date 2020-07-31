# Tools

* "fish plot": [TimeScape](http://bioconductor.org/packages/release/bioc/html/timescape.html)
* clonality inference: [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Home) + [CITUP](https://github.com/amcpherson/citup)
* clonality inference: PhyloWGS, DPClust, CliP?

## single cell CNV analysis:

* pipeline: [Genome-wide copy number analysis of single cells (2012)](https://www.nature.com/articles/nprot.2012.039) and
[Optimizing sparse sequencing of single cells for highly multiplex copy number profiling (2015)](https://genome.cshlp.org/content/25/5/714)
* DNAcopy: used for circular binary segmentation of CNV bins, then `MergeLevels` can be used to join similar segments [1]
* dbscan used to filter noise [1]
* 'cluster' package: to find optimal number of clusters [1]
* hierarchical clustering: R `stats::hclust()` [1] (ward.D2 method)
* measure internal cluster consistency using Pearson and Spearman correlations

## Mutations analysis:

* 2-dimensional clustering: MDS, r function `cmdscale(x, eig = TRUE, k = 2)` where x is a matrix with cells in the columns and mutations in the rows

## comparison of phylogenetic trees:

* Tanglegrams, `dendextend` R package


1. [Multiclonal Invasion in Breast Tumors Identified by Topographic Single Cell Sequencing, Navin 2018](https://www.sciencedirect.com/science/article/pii/S0092867417314496)

# Notes from articles

## Single-cell DNA sequencing reveals a latedissemination model in metastatic colorectal cancer, M.L. Leung et al. (2017)

* use PyClone to normalize the VAFs by copy number events and calculated clonal frequencies to check if VAF differences between 
primary tumor and metastases are caused by selection or CNV event
* identify subpopulations of cells: perform MDS scaling of mutation data and cluster cells, plot dim 1 vs dim 2 (`cmdscale(x, eig = TRUE, k = 2)`)
* continue with hierarchical clustering using heatmaps (`heatmap.2()` from qplots (CRAN))
* do MDS and hierarhical clustering also for CNV data

    >Pairwise Euclidean distances were calculated from the single-cell
    >copy number data matrix (log2[ratio + 0.1]) and then used for hierarchical
    >clustering using ward-linkage in R using the heatmap.3
    >function from the “gplots” package available on CRAN
    
* use SCITE (Jahn et al. 2016) to compute phylogenetic mutation trees, use Cytoscape to redraw it
* then they integrate CNV and mutatio trees (how?)
* VCF annotation: ANNOVAR (Wang et al. 2010)
