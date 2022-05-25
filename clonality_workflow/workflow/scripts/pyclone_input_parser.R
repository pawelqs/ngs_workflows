### CURRENTLY NOT USED
### This script generates all files (both tsv and yaml) needed by PyClone
### Unfortunatelly, I can not create sample-wise tsv and patient-wise yaml in the same Snakemake rule
### Therefore the secendo script was preapred, which generates only sample-tsv files
### and config.yaml for patient must be created in another rule

suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(tidyverse))

###############################
###### Parse args
###############################

option_list = list(
  make_option(c("--vcf_file"), action="store",
              help="Multisample VCF file"),
  make_option(c("--cnv_files"), action="store",
              help="Comma-separated list of files with CNV segments for each file"),
  make_option(c("--sample_names"), action="store",
              help="Comma-separated list of samples"),
  make_option(c("--cnv_type"), action="store",
              help="CNV detection algorithm. Supprted: facets"),
  make_option(c("--working_dir"), action="store",
              help="Working directory")  ,
  make_option(c("--sex"), action="store",
              help="Patient sex")  ,
  make_option(c("--genome"), action="store",
              help="Genome version")  
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

vcf_file <- opt$vcf
cnv_files <- str_split(opt$cnv_files, ",")[[1]]
working_dir <- opt$working_dir
samples <- str_split(opt$sample_names, ",")[[1]]
cnv_type <- opt$cnv_type
sex <- opt$sex
genome <- opt$genome

# samples <- c("G30_P1", "G30_L1")
# vcf_file <- "results/filtered_vcf/G30.vcf"
# cnv_files <- str_c("results/facets/", samples, ".csv")
# working_dir <- "results/pyclone_facets/G30"
# sex <- "female"
# genome <- "hg38"
# cnv_type <- "facets"

names(samples) <- samples
names(cnv_files) <- samples

#####################################
############# Prepare output
#####################################

tsv_dir <- str_c(working_dir, "/input") %>%
  str_replace("//", "/")
if (!dir.exists(tsv_dir)) {
  dir.create(tsv_dir, recursive = TRUE)
}
out_config <- str_c(working_dir, "/config.yml") %>%
  str_replace("//", "/")
out_tsv <- str_c(tsv_dir, "/", samples, ".tsv")
names(out_tsv) <- samples

#####################################
############# Prepare tsv files
#####################################

chromosomes <- str_c("chr", c(1:22, "X", "Y"))

parse_facets_csv <- function(cnv_file) {
  read_csv(cnv_file) %>%
    transmute(
      seqnames = str_c("chr", chrom) %>%
        str_replace("chr23", "chrX") %>%
        str_replace("chr24", "chrY"),
      start = start,
      end = end,
      minor_cn = lcn.em,
      major_cn = tcn.em - lcn.em
    )
}

create_cnv_granges <- function(cnv, genome) {
  seq_info <- plyranges::genome_info(genome) %>%
    seqinfo() %>% 
    as.data.frame() %>%
    rownames_to_column("chr") %>%
    filter(chr %in% chromosomes)

  as_granges(cnv) %>%
    set_genome_info( 
      genome = "hg38", 
      seqnames = seq_info$chr,
      seqlengths = seq_info$seqlengths,
      is_circular = seq_info$isCircular
    )
}

add_normal_regions <- function(cnv) {
  normal_regions <- cnv %>%
    compute_coverage() %>%
    filter(score == 0) %>%
    select(-score)
  
  c(cnv, normal_regions) %>%
    sort() %>%
    mutate(
      sex,
      chr = as.character(seqnames),
      normal_cn = if_else(chr == "chrY" & sex == "male", 1, 2),
      minor_cn = replace_na(minor_cn, 0),
      major_cn = if_else(is.na(major_cn), normal_cn, major_cn)
    )
}

prepare_mutation_data <- function(sample, tidy_vcf) {
  inner_join(
      filter(tidy_vcf$fix, FILTER == "PASS"),
      filter(tidy_vcf$gt, Indiv == sample, !is.na(gt_AD))
    ) %>%
    transmute(
      seqnames = CHROM,
      start = POS,
      end = POS,
      mutation_id = str_c(seqnames, POS, sep = ":"),
      gt_AD,
      n_variants = str_split(gt_AD, ",") %>% map_int(length)
    ) %>%
    filter(n_variants == 2) %>%
    separate(gt_AD, into = c("ref_counts", "var_counts"))
}

join_mutations_with_cnv <- function(muts, cnvs) {
  muts %>%
    as_granges() %>%
    join_overlap_left(cnvs) %>%
    as.data.frame() %>%
    select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn) %>%
    drop_na() %>%
    filter(major_cn > 0)
}


tidy_vcf <- vcf_file %>%
  read.vcfR() %>%
  vcfR2tidy(format_fields = "AD")

cnvs_parsed <- map(cnv_files, parse_facets_csv)
cnvs_ranges <- map(cnvs_parsed, create_cnv_granges, genome)
cnvs_ranges2 <- map(cnvs_ranges, add_normal_regions)

muts <- map(samples, prepare_mutation_data, tidy_vcf)
mut_data <- map2(muts, cnvs_ranges2, join_mutations_with_cnv)

walk2(mut_data, out_tsv, write_tsv)

########################################
######### Prepare config files
########################################

get_sample_purity <- function(sample, cnv_type) {
  if (cnv_type == "facets") {
    cnv_files[sample] %>%
      read_csv() %>%
      pull(Purity) %>%
      unique()
  } else 1
}
  
samples_data <- map(samples, function(sample) {
  list(
    mutations_file = out_tsv[sample],
    tumour_content = list(
      value = get_sample_purity(sample, cnv_type)
    ),
    error_rate = 0.001
  )
})

config_yml <- list(
  working_dir = working_dir,
  trace_dir = "trace/pyclone_binomial/all",
  density = "pyclone_binomial",
  num_iters = 10000L,
  base_measure_params = list(
    alpha = 1L,
    beta = 1L
  ),
  concentration = list(
    value = 1.0
  ),
  prior = list(
    shape = 1.0,
    rate = 0.001
  ),
  samples = samples_data
)
print(config_yml)
print(out_config)
write_yaml(config_yml, out_config)
