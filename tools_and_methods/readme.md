# Tools

* "fish plot": [TimeScape](http://bioconductor.org/packages/release/bioc/html/timescape.html)
* clonality interence: [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Home) + [CITUP](https://github.com/amcpherson/citup)

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
