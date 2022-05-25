library(sciClone, quietly = TRUE)
library(vcfR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(tidyverse, quietly = TRUE)

###############################
###### Parse args
###############################

parse_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- args %>%
    str_c(collapse = " ") %>%
    str_split("--") %>%
    unlist() %>%
    .[-1] %>%
    str_trim() %>%
    str_split(" ")
  opt <- map(args, ~.x[-1])
  names(opt) <- map(args, 1)
  return(opt)
}

usage <- "Rscript.exe run_sciclone.R 
  --vcf G30.vcf 
  --cnv_files G30_P1.csv G30_L1.csv 
  --sample_names G30_P1 G30_L1 
  --cnv_type FACETS 
  --prefix G30 
  --target_regions targets.bed
  --outdir results/sciclone\n
"

opt <- parse_arguments()
required_fields <- c("vcf", "cnv_files", "sample_names", "cnv_type", "prefix", "outdir")
opt_null <- map_lgl(required_fields, ~is.null(opt[[.x]]))
if(any(opt_null)) {
  cat("\n--------------------USAGE--------------------\n\n")
  cat(usage)
  stop()
}

vcf_file <- opt$vcf
cnv_files <- opt$cnv_files
samples <- opt$sample_names
cnv_type <- opt$cnv_type
patient <- opt$prefix
taget_regions_file <- opt$target_regions
out_dir <- opt$outdir

print(opt)

# patient <- "G30"
# samples <- c("G30_P1", "G30_L1")
# vcf_file <- str_c("results/filtered_vcf/", patient, ".vcf")
# cnv_files <- str_c("results/facets/", samples, ".csv")
# cnv_type <- "FACETS"
# taget_regions_file <- "resources/agilent_sureselect_v6/Agilent_SureSelect_v6r2_S07604514_Covered_hg38.cnvkit.bed"
# out_dir <- "temp"

names(samples) <- samples
names(cnv_files) <- samples


###############################
###### prepare output
###############################

out_dir <- str_replace(out_dir, "/$", "")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- str_c(out_dir, "/", patient, ".sciclone.tsv")
out_1d <- str_c(out_dir, "/", patient, ".sciclone_1d.pdf")
out_2d <- str_c(out_dir, "/", patient, ".sciclone_2d.pdf")
out_3d <- str_c(out_dir, "/", patient, ".sciclone_3d.pdf")
out_stats <- str_c(out_dir, "/", patient, ".sciclone.stats")


###############################
###### Prepare input data
###############################

vcf <- read.vcfR(vcf_file)
tidy_vcf <- vcfR2tidy(vcf, format_fields = c("AD", "DP", "AF"))

# 5 column, tab delimited: chr, pos, ref_reads, var_reads, vaf
get_mut_data <- function(sample) {
  inner_join(
    tidy_vcf$fix, 
    filter(tidy_vcf$gt, Indiv == sample)
  ) %>%
    select(chr = CHROM, st = POS, gt_AD, vaf = gt_AF) %>%
    separate(gt_AD, into = c("ref_reads", "var_reads")) %>%
    mutate(
      chr,
      ref_reads = parse_integer(ref_reads),
      var_reads = parse_integer(var_reads),
      vaf = parse_double(vaf)
    ) %>%
    as.data.frame()
}
mut_data <- map(samples, get_mut_data)

# 4 columns - chr, start, stop, segment_mean
get_cnv_data <- function(file, file_type) {
  if (file_type == "FACETS") {
    # segment median used instead of seg mean!!!!
    read_csv(file) %>%
      transmute(
        chr = str_c("chr", chrom) %>%
          str_replace("chr23", "chrX") %>%
          str_replace("chr24", "chrY"),
        start,
        stop = end,
        segment_median = cnlr.median
      ) %>%
      as.data.frame()
  }
}
cnv_data <- map(cnv_files, get_cnv_data, file_type = cnv_type)

# 3 column BED
get_LOHs <- function(file, file_type) {
  if (file_type == "FACETS") {
    read_csv(file, col_types = cols(
        chrom = col_character(),
        start = col_integer(),
        end = col_integer(),
        tcn.em = col_integer(),
        lcn.em = col_integer()
        )
      ) %>%
      filter(lcn.em == 0) %>%
      mutate(
        seqnames = str_c("chr", chrom) %>%
          str_replace("chr23", "chrX") %>%
          str_replace("chr24", "chrY")
      ) %>%
      # select(chrom, start, end, tcn.em, lcn.em) %>%
      select(seqnames, start, end) %>%
      plyranges::as_granges()
  }
}
LOHs <- map(cnv_files, get_LOHs, file_type = "FACETS")

exclude_ranges <- LOHs %>%
  GRangesList() %>%
  unlist() %>%
  plyranges::reduce_ranges()

target_ranges <- read_tsv(
    taget_regions_file, 
    skip = 2, 
    col_names = c("seqnames", "start", "end", "."),
    col_types = cols_only(
      seqnames = col_character(), 
      start = col_integer(), 
      end = col_integer()
    )
  ) %>%
  plyranges::as_granges()

target_length <- as_tibble(target_ranges) %>%
  pull(width) %>%
  sum

exclude_length <- plyranges::intersect_ranges(target_ranges, exclude_ranges) %>%
  as_tibble() %>%
  pull(width) %>%
  sum

exclude_reg <- as.data.frame(exclude_ranges) %>%
  select(seqnames, start, end) %>%
  mutate(seqnames = as.character(seqnames))

msg <- sprintf("%.2f\tproportion of the targeted sequences excluded as LOH regions", 
            exclude_length / target_length)
write_file(msg, file = out_stats)

###############################
###### Run sciClone
###############################

safe_sciClone <- safely(sciClone)

sc <- safe_sciClone(
  vafs = mut_data,
  copyNumberCalls = cnv_data,
  sampleNames = samples,
  regionsToExclude = exclude_reg
)

if (!is.null(sc$result)) {
  sc <- sc$result
  writeClusterTable(sc, out_file)
  sc.plot1d(sc, out_1d)
  if(length(samples) >=2) {
    sc.plot2d(sc, out_2d)
  }
  if(length(samples) >=3) {
    sc.plot3d(sc, sc@sampleNames, size=700, outputFile = out_3d)
  }
} else {
  write_file("ERROR", out_file)
}
