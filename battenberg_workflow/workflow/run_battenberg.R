library(Battenberg)

opt <- list(
    nb = snakemake@input[["nb"]],
    tb = snakemake@input[["tb"]],
    normalname = snakemake@params[["normalname"]],
    tumoruname = snakemake@params[["tumourname"]],
    sex = snakemake@params[["sex"]],
    output = snakemake@output[[1"output"]],
    skip_allelecount = snakemake@params[["skip_allelecount"]],
    skip_preprocessing = snakemake@params[["skip_preprocessing"]],
    skip_phasing = snakemake@params[["skip_phasing"]],
    cpu = snakemake@params[["cpu"]],
    bp = snakemake@params[["bp"]],
    outfile = snakemake@output[[1]]
)

# opt$tb <- strsplit(opt$tb, ",")[[1]]
# opt$tumourname <- strsplit(opt$tumourname, ",")[[1]]

print(opt)

# TUMOURNAME = opt$tumourname
# NORMALNAME = opt$normalname
# NORMALBAM = opt$nb
# TUMOURBAM = opt$tb
# IS.MALE = opt$sex=="male" | opt$sex=="Male"
# RUN_DIR = opt$output
# SKIP_ALLELECOUNTING = opt$skip_allelecount
# SKIP_PREPROCESSING = opt$skip_preprocessing
# SKIP_PHASING = opt$skip_phasing
# NTHREADS = opt$cpu
# PRIOR_BREAKPOINTS_FILE = opt$bp

# analysis = "paired"

# JAVAJRE = "java"
# ALLELECOUNTER = "alleleCounter"
# IMPUTE_EXE = "impute2"

# GENOMEBUILD = "hg38"
# USEBEAGLE = T

# BEAGLE_BASEDIR = "/home/pkus/resources/battenberg_hg38_new/chrX"
  
# IMPUTEINFOFILE = file.path(BEAGLE_BASEDIR, "impute_info.txt")
# G1000ALLELESPREFIX = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_")
# G1000LOCIPREFIX = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_")
# GCCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "GC_correction_hg38/1000G_GC_")
# REPLICCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "RT_correction_hg38/1000G_RT_")
# PROBLEMLOCI = file.path(BEAGLE_BASEDIR, "probloci/probloci.txt.gz")
	
# BEAGLEREF.template = file.path(BEAGLE_BASEDIR, "beagle/CHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz")
# BEAGLEPLINK.template = file.path(BEAGLE_BASEDIR, "beagle/plink.CHROMNAME.GRCh38.map")
# BEAGLEJAR = file.path(BEAGLE_BASEDIR, "beagle/beagle.08Feb22.fa4.jar")

# # CHROM_COORD_FILE = "/homes/sdentro/repo/battenberg/chromosome_coordinates_hg38.txt"

# PLATFORM_GAMMA = 1
# PHASING_GAMMA = 1
# SEGMENTATION_GAMMA = 10
# SEGMENTATIIN_KMIN = 3
# PHASING_KMIN = 1
# CLONALITY_DIST_METRIC = 0
# ASCAT_DIST_METRIC = 1
# MIN_PLOIDY = 1.6
# MAX_PLOIDY = 4.8
# MIN_RHO = 0.1
# MIN_GOODNESS_OF_FIT = 0.63
# BALANCED_THRESHOLD = 0.51
# MIN_NORMAL_DEPTH = 10
# MIN_BASE_QUAL = 20
# MIN_MAP_QUAL = 35
# CALC_SEG_BAF_OPTION = 1


# Change to work directory and load the chromosome information
# setwd(RUN_DIR)

# battenberg(
#   analysis=analysis,
#   tumourname=TUMOURNAME, 
#   normalname=NORMALNAME, 
#   tumour_data_file=TUMOURBAM, 
#   normal_data_file=NORMALBAM, 
#   ismale=IS.MALE, 
#   imputeinfofile=IMPUTEINFOFILE, 
#   g1000prefix=G1000LOCIPREFIX, 
#   g1000allelesprefix=G1000ALLELESPREFIX, 
#   gccorrectprefix=GCCORRECTPREFIX, 
#   repliccorrectprefix=REPLICCORRECTPREFIX, 
#   problemloci=PROBLEMLOCI, 
#   data_type="wgs",
#   impute_exe=IMPUTE_EXE,
#   allelecounter_exe=ALLELECOUNTER,
#   usebeagle=USEBEAGLE,
# 	beaglejar=BEAGLEJAR,
# 	beagleref=BEAGLEREF.template,
# 	beagleplink=BEAGLEPLINK.template,
# 	beaglemaxmem=10,
# 	beaglenthreads=1,
# 	beaglewindow=40,
# 	beagleoverlap=4,
# 	javajre=JAVAJRE,
#   nthreads=NTHREADS,
#   platform_gamma=PLATFORM_GAMMA,
#   phasing_gamma=PHASING_GAMMA,
#   segmentation_gamma=SEGMENTATION_GAMMA,
#   segmentation_kmin=SEGMENTATIIN_KMIN,
#   phasing_kmin=PHASING_KMIN,
#   clonality_dist_metric=CLONALITY_DIST_METRIC,
#   ascat_dist_metric=ASCAT_DIST_METRIC,
#   min_ploidy=MIN_PLOIDY,
#   max_ploidy=MAX_PLOIDY,
#   min_rho=MIN_RHO,
#   min_goodness=MIN_GOODNESS_OF_FIT,
#   uninformative_BAF_threshold=BALANCED_THRESHOLD,
#   min_normal_depth=MIN_NORMAL_DEPTH,
#   min_base_qual=MIN_BASE_QUAL,
#   min_map_qual=MIN_MAP_QUAL,
#   calc_seg_baf_option=CALC_SEG_BAF_OPTION,
#   skip_allele_counting=SKIP_ALLELECOUNTING,
#   skip_preprocessing=SKIP_PREPROCESSING,
#   skip_phasing=SKIP_PHASING,
#   prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
#   GENOMEBUILD=GENOMEBUILD
# 	# chrom_coord_file=CHROM_COORD_FILE
# )

outfile <- file(outfile)
writeLines(c("Hello","World"), outfile)
close(outfile)