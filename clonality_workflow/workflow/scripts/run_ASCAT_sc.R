library(ASCAT.sc)


bams <- snakemake@input[["bams"]]
output <- snakemake@output[["rda"]] ## outputs here
genome_fa <- snakemake@params[["genome_fa"]]
workdir <- snakemake@params[["workdir"]] ## outputs here
genome <- snakemake@params[["genome"]]
projectname <- snakemake@params[["projectname"]]
multipcf <- snakemake@params[["multipcf"]]
sex <- snakemake@params[["sex"]]
chrstring_bam <- snakemake@params[["chrstring_bam"]]
threads <- snakemake@threads

dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

res <- run_sc_sequencing(
    tumour_bams = bams,
    allchr = paste0(chrstring_bam, c(1:22, "X", "Y")), ## might need a "chr" instead of "" if hg38
    sex = rep(sex, length(bams)),
    chrstring_bam = chrstring_bam,          ## might need a "chr" instead of "" if hg38
    purs = seq(0.1, 1, 0.01) ,              ## purity grid values
    ploidies = seq(1.7,5, 0.01),            ## average ploidy grid values
    maxtumourpsi = 5,                       ## maximum tumour ploidy
    binsize = 500000,                       ## bin size - reduce if enough sequencing reads (look at dpb (depth per bin) value in plots can go down to 100bp or even lower)
    build = genome,                         ## either hg19 or hg38 so far
    MC.CORES = threads,                     ## number of cores available
    outdir = workdir,
    projectname = projectname,
    multipcf = multipcf                     ## lgl use multipcf for multi-track segmentation if multi-sample sequencing
)

saveRDS(res, file = output)
