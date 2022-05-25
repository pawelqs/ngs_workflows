rm(list = ls())
library(jsonlite)
library(tidyverse)
library(optparse)


###############################
###### Parse args & prepare names
###############################

option_list = list(
  make_option("--summ_json_file", action="store", help="Gunzipped json file"),
  make_option("--muts_json_file", action="store", help="Gunzipped json file"),
  make_option("--mutass_zip_file", action="store", help="ZIP file"),
  make_option("--out_dir", action="store", help="Output dir"),
  make_option("--prefix", action="store", help="Prefix for output files eg. patient id")
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

summ_json_file <- opt$summ_json_file
muts_json_file <- opt$muts_json_file
mutass_zip_file <- opt$mutass_zip_file
out_dir <- opt$out_dir
prefix <- opt$prefix


# prefix <- "G30_facets"
# out_dir <- "."
# patient_dir <- "../data_raw/pwgs_titan_2021_02_01/G30"
# summ_json_file <- str_c(patient_dir, ".summ.json")
# muts_json_file <- str_c(patient_dir, ".muts.json")
# mutass_zip_file <- str_c(patient_dir, ".mutass.zip")

# prepare output files
temp_dir <- str_c(prefix, "_temp")
rds_file <- str_c(out_dir, "/", prefix, ".Rds") %>% str_replace("//", "/")
tree_file <- str_c(out_dir, "/", prefix, ".tree.tsv") %>% str_replace("//", "/")
clusters_file <- str_c(out_dir, "/", prefix, ".clusters.tsv") %>% str_replace("//", "/")
mutations_file <- str_c(out_dir, "/", prefix, ".mutations.tsv") %>% str_replace("//", "/")


###############################
###### Extract best
###############################

  results <- list()
  
  ####################### Find the best tree
  # best tree according to max tree_density - like in PhyloWGS/witness
  summ <- summ_json_file %>% read_json()
  tree_densities <- tibble(
    tree = names(summ$tree_densities),
    density = unlist(summ$tree_densities)
  )
  best_tree_id <- tree_densities %>%
    arrange(desc(density)) %>%
    slice(1) %>%
    pull(tree)
  
  ## Theoretically also good way, but gives a different result is:
  # llh <- map_dbl(summ$trees, "llh")
  # best_tree_id <- which.max(llh)
  
  ####################### Get tree details
  results$clustering_index <- summ$trees[[best_tree_id]]$clustering_index
  results$branching_index <- summ$trees[[best_tree_id]]$branching_index
  results$linearity_index <- summ$trees[[best_tree_id]]$linearity_index
  results$llh <- summ$trees[[best_tree_id]]$llh
  
  ####################### Get tree details
  results$node_childs <- summ$trees[[best_tree_id]]$structure %>% enframe
  
  ####################### Get clusters' details
  results$clusters_info <- summ$trees[[best_tree_id]]$populations %>%
    imap(function(population, cluster_id) {
      tibble(
        cluster = cluster_id,
        num_ssms = population$num_ssms,
        cellularities = list(population$cellular_prevalence)
      )
    }) %>%
    bind_rows()
  
  ####################### Get mutations' details
  muts <- muts_json_file %>% read_json()
  mutations <- muts$ssms %>%
    imap(function(x, name) {
      res <- list(
        pwgs_mut_id = name,
        id = x$name,
        total_reads = list(x$total_reads),
        ref_reads = list(x$ref_reads)
      )
    }) %>%
    bind_rows() %>%
    rename(mutation_id = id)
  
  mutass <- unzip(
    zipfile = mutass_zip_file,
    files = str_c(best_tree_id, ".json"),
    exdir = temp_dir
    ) %>% read_json()
  clusters <- mutass$mut_assignments %>%
    imap(function(m, cluster) {
      tibble(
        cluster = cluster,
        pwgs_mut_id = unlist(m)
      )
    }) %>%
    bind_rows()
  
  results$mutations <- inner_join(clusters, mutations, by = "pwgs_mut_id")

###############################
###### Save results
###############################

saveRDS(results, file = rds_file)

stringify_list <- function(x) {
  x %>%
    map(unlist) %>%
    map_chr(str_c, collapse = ",")
}

results$node_childs %>%
  mutate_if(is_list, stringify_list) %>%
  write_tsv(file = tree_file)
  
results$clusters_info %>%
  mutate_if(is_list, stringify_list) %>%
  write_tsv(file = clusters_file)

results$mutations %>%
  mutate_if(is_list, stringify_list) %>%
  write_tsv(file = mutations_file)

unlink(temp_dir, recursive = TRUE)