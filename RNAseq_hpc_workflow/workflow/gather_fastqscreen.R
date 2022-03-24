library(tidyverse)

# Read files
# files <- list.files("../QC/trimmed_screen/", pattern = ".txt", full.names = TRUE)
files <- snakemake@input %>%
  unlist() %>%
  sort()
names(files) <- files %>%
  str_replace("qc/trimmed/", "") %>%
  str_replace("_val_1_screen.txt", "") %>%
  str_replace("_val_2_screen.txt", "")
print(files)

dt <- map(files, read_tsv, skip = 1) %>%
  bind_rows(.id = "file")




cols <- tibble(
  levs = c("%One_hit_one_genome", "%Multiple_hits_one_genome", "%One_hit_multiple_genomes",
           "%Multiple_hits_multiple_genomes", "%Unmapped"),
  colors = rev(c("royalblue4", "royalblue", "red", "red4", "gray"))
)

w <- 15
rows <- ceiling(length(files)/6)
h <- 1 + rows*(7/6)

p <- dt %>%
  select(-starts_with("#"), -Multiple_hits_multiple_genomes) %>%
  gather(key = "stat", value = "values", `%Unmapped`:`%Multiple_hits_multiple_genomes`) %>%
  filter(!str_detect(Genome, "%")) %>%
  mutate(stat = parse_factor(stat, levels = rev(cols$levs))) %>%
  ggplot(aes(Genome, values, fill = stat)) +
    geom_bar(stat = "identity") +
    facet_wrap(~file, ncol = 6) +
    coord_flip() +
    scale_fill_manual(values = cols$colors, breaks = cols$levs) +
    theme(legend.position = "top") +
    guides(fill = guide_legend(ncol = 3))
ggsave(filename = snakemake@output[[1]], plot = p, width = w, height = h)
