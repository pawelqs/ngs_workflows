suppressPackageStartupMessages(library(clonalityParsers))

prepare_pycloneVI_input(
  vcf_file = snakemake@input[["vcf"]],
  cnv_files = snakemake@input[["cnv"]],
  sample_ids = snakemake@params[["sample_ids"]], 
  sex = snakemake@params[["sex"]], 
  genome_build = snakemake@params[["genome_build"]],
  snv_algorithm = snakemake@params[["snv_algorithm"]],
  cnv_algorithm = snakemake@params[["cnv_algorithm"]],
  filename = snakemake@output[[1]]
)
