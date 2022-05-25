suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(tidyverse))

###############################
###### Parse args
###############################

option_list = list(
  make_option("--vcf_file", action="store", help="VCF file"),
  make_option("--cnv_file", action="store", help="CNV file"),
  make_option("--out_file", action="store", help="TSV output file"),
  make_option("--sample_name", action="store", help="Name of the sample in VCF file"),
  make_option("--cnv_type", action="store", help="CNV detection algorithm. Supprted: facets"),
  make_option("--sex", action="store", help="Patient sex")  ,
  make_option("--genome", action="store", help="Genome version")  
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

vcf_file <- opt$vcf
cnv_file <- opt$cnv_file
tsv_file <- opt$out_file
sample <- opt$sample_name
cnv_type <- opt$cnv_type
sex <- opt$sex
genome <- opt$genome

# vcf_file <- "results/filtered_vcf/AMLRO-9.vcf"
# cnv_file <- "results/facets/AMLRO-9_Dx.csv"
# tsv_file <- "results/pyclone_facets/AMLRO-9/input/AMLRO-9_Dx.tsv"
# sample <- "AMLRO-9_Dx"
# cnv_type <- "facets"
# sex <- "male"
# genome <- "hg38"

#####################################
############# Prepare output
#####################################

tsv_dir <- dirname(tsv_file)
if (!dir.exists(tsv_dir)) {
  dir.create(tsv_dir, recursive = TRUE)
}

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
    ) %>%
    drop_na()
}

parse_titan_tsv <- function(cnv_file) {
  read_tsv(cnv_file) %>%
    transmute(
      seqnames = str_c("chr", Chromosome) %>%
        str_replace("chr23", "chrX") %>%
        str_replace("chr24", "chrY"),
      start = `Start_Position(bp)`,
      end = `End_Position(bp)`,
      minor_cn = MinorCN,
      major_cn = MajorCN
    ) %>%
    drop_na()
}

create_cnv_granges <- function(cnv, genomeBuild) {
  seq_info <- plyranges::genome_info(genomeBuild) %>%
    seqinfo() %>% 
    as.data.frame() %>%
    rownames_to_column("chr") %>%
    filter(chr %in% chromosomes)
  
  as_granges(cnv) %>%
    set_genome_info( 
      genome = genomeBuild, 
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
      mutation_id = str_c(seqnames, POS, sep = "_"),
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

cnvs_parsed <- if (cnv_type == "facets") {
    parse_facets_csv(cnv_file)
  } else if (cnv_type == "titan") {
    parse_titan_tsv(cnv_file)
  }
cnvs_ranges <- create_cnv_granges(cnvs_parsed, genome)
cnvs_ranges2 <- add_normal_regions(cnvs_ranges)

muts <- prepare_mutation_data(sample, tidy_vcf)
mut_data <- join_mutations_with_cnv(muts, cnvs_ranges2)

write_tsv(mut_data, tsv_file)
