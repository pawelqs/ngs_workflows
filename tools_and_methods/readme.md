# Tools

* VCF annotation: ANNOVAR (Wang et al. 2010)

## Reconstruction of subclonal structure:

* [SCHISM](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004416)
* [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Home) + [CITUP](https://github.com/amcpherson/citup)
* PhyloWGS
* DPClust
* CliP
* "fish plot": [TimeScape](http://bioconductor.org/packages/release/bioc/html/timescape.html)

## comparison of phylogenetic trees:


# Tools for single cells

[list of over 100 tools](https://github.com/seandavi/awesome-single-cell)

## Identification of SNV:

*[Identification of somatic mutations in single cell DNA-seq using a spatial model of allelic imbalance](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6715686/)
	* Monovar (makes corrections for for allelic dropout)
	* SCcaller - they developped SCMDA
	* SCAN-SNV - for MDA amplification results <- their method, do not take into account CNV

## Identification of CNV:

* Varbin algorithm:

	> The sequencing data was processed following the ‘variable binning’ pipeline (Baslan et al., 2012; Baslan et al., 2015). Briefly, reads were aligned to the human genome HG19 using Bowtie2 and counted in variable bins at a genomic resolution of 220kb. Unique normalized read counts were segmented using the circular binary segmentation (CBS) method from R Bioconductor ‘DNAcopy’ package (Shah et al., 2006) followed by MergeLevels to join adjacent segments with non-significant differences in segmented ratios. The parameters used for CBS segmentation were alpha=0.0001 and undo.prune=0.05. Default parameters were used for MergeLevels, which removed erroneous chromosome breakpoints. Data was filtering with more than 100 break points or identified as noise with the R package for Density-based spatial clustering of applications (Dudik et al., 2015; Martin Ester, 1996; Yue et al., 2004) with noise ‘dbscan’ (v1.1-1) (Piekenbrock, 2017). We used this package to determine technical noise within the copy number profiles. We examined the plots to find the elbow and recorded this value for selecting the eps to filter data using the ‘dbscan’ package. Using this number, dbscan determined which single cell samples exhibited to much technical noise and filtered approximately 20% of the total datasets for each patient.

* [Assessing the performance of methods for copy number aberration detection from single-cell DNA sequencing data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7377518/)
	* HMMcopy
	* CopyNumber
	* Ginkgo (winner)

## CNV analysis:

* [methods comparison by Navin](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7377518/)
* pipeline: [Genome-wide copy number analysis of single cells (2012)](https://www.nature.com/articles/nprot.2012.039) and
[Optimizing sparse sequencing of single cells for highly multiplex copy number profiling (2015)](https://genome.cshlp.org/content/25/5/714)
* DNAcopy: used for circular binary segmentation of CNV bins, then `MergeLevels` can be used to join similar segments [1]
* dbscan used to filter noise [1]
* 'cluster' package: to find optimal number of clusters [1]
* hierarchical clustering: R `stats::hclust()` [1] (ward.D2 method)
* measure internal cluster consistency using Pearson and Spearman correlations

## mutations analysis:

* To test whether the genotype matrix violates the infinite-sites assumption, we ran the four-gamete test. The four-gamete theorem states that an m × n binary matrix M has 
an undirected perfect phylogeny if and only if no pair of columns contain all four binary pairs (0,0; 0,1; 1,0; and 1,1), where m represents the number of taxa (leaves of the tree) 
and n represents genomic sites [47].

* 2-dimensional clustering in the MDS space (r function `cmdscale(x, eig = TRUE, k = 2)` where x is a matrix with cells in the columns and mutations in the rows)
* continue with hierarchical clustering using heatmaps (`heatmap.2()` from qplots (CRAN))

**Phylogenetic Inference:**
* use SCITE (Jahn et al. 2016) to compute phylogenetic mutation trees, use Cytoscape to redraw it [1]
* OncoNEM: (Ross et al. 2016)
* Single Cell Genotyper: (Roth et al. 2016)
* SiFit (Zafar et al. 2017) - uses finite site model

**comparison of tree**
* False negative (FN) distance : This counts the edges in Tt that induce bipartitions that are not present in C(Ti). This distance is normalized by dividing by the total number of bipartitions in Tt,
* False positive (FP) distance : This counts the edges in Ti that induce bipartitions that are not present in C(Tt). This distance is normalized by dividing by the total number of bipartitions in Ti
* Robinson–Foulds (RF) distance : The Robinson–Foulds distance is the average of the FP and FN distances. This is the most common error metric `dist.topo()` from `ape` package
* `dendextend` R package
	* tanglegram()
	* dend_diff()
	* entanglement()
	* cor_bakers_gamma() - Baker’s Gamma Index

**Shah Lab**

* TitanCNA
* HMMcopy


**R packages**

* tidyverse
* GenomicRanges, plyranges
* DESeq2. edgeR
* fgsea, tmod

# Wiki

## Methods of single cell genome amplification:

* MDA
* MALBAC
* DOP-PCR
* TnBC
-> https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181163

# Articles

1. [Multiclonal Invasion in Breast Tumors Identified by Topographic Single Cell Sequencing, Navin 2018](https://www.sciencedirect.com/science/article/pii/S0092867417314496)

2. Single-cell DNA sequencing reveals a latedissemination model in metastatic colorectal cancer, M.L. Leung et al. (2017)

	* use PyClone to normalize the VAFs by copy number events and calculated clonal frequencies to check if VAF differences between 
	primary tumor and metastases are caused by selection or CNV event
	* do MDS and hierarhical clustering also for CNV data
		>Pairwise Euclidean distances were calculated from the single-cell
		>copy number data matrix (log2[ratio + 0.1]) and then used for hierarchical
		>clustering using ward-linkage in R using the heatmap.3
		>function from the “gplots” package available on CRAN
	* then they integrate CNV and mutatio trees (how?)
