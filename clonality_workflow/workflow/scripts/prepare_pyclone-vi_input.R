# install.packages("tidyverse")
# install.packages("remotes")
# install.packages("optparse")

remotes::install_github("https://github.com/pawel125/clonalityParsers.git")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(clonalityParsers))
suppressPackageStartupMessages(library(stringr))

option_list <- list(
  make_option(c("-s", "--vcf_file"), dest = "vcf_file",
              help="Path to VCF file"),
  make_option(c("-c", "--cnv_files"), dest="cnv_files",
              help="Path to CNV files (multiple values migh be comma-separated"),
  make_option(c("-i", "--sample_ids"), dest = "sample_ids",
              help="Sample IDs"),
  make_option("--sex", dest = "sex",
              help = "male/female"),
  make_option("--genome_build", dest = "genome_build",
              help="eg. hg38, passed to GenomeInfoDb::seqinfo()"),
  make_option("--snv_algorithm", dest = "snv_algorithm",
              help="algorithm used for SNVs detection, will be recognized if NULL, default: NULL"),
  make_option("--cnv_algorithm", dest = "cnv_algorithm",
              help="algorithm used for CNVs detection, will be recognized if NULL, default: NULL"),
  make_option("--filename", dest = "filename",
              help="File name to create on disk, default: NULL, no save")
)

opt <- parse_args(OptionParser(option_list=option_list))
opt$cnv_files <- str_split(opt$cnv_files, pattern = ",")[[1]]
opt$sample_ids <- str_split(opt$sample_ids, pattern = ",")[[1]]
print(opt)

prepare_pycloneVI_input(
  vcf_file = opt$vcf_file,
  cnv_files = opt$cnv_files,
  sample_ids = opt$sample_ids, 
  sex = opt$sex, 
  genome_build = opt$genome_build,
  snv_algorithm = opt$snv_algorithm,
  cnv_algorithm = opt$cnv_algorithm,
  filename =  = opt$filename
)

