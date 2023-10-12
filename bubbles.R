library(tidyverse)
library(ComplexHeatmap)

theme_set(theme_bw())

nodes <- read_delim("gfa/graph_nodes.csv", delim = " ", col_names = c("node", "type"))

cpgs <- read_delim("gfa/cpg_index_xaf.csv.gz", delim = " ", col_names = c("node", "pos"))

load_alt_cpg <- function(path, pattern, nodes, cpgs) {
  cpg <- bind_rows(lapply(
    list.files("data/bubbles/", pattern = pattern, full.names = T),
    function(x) {
      read_delim(x,
        col_names = c("node", "pos", "strand", "mc", "cov", "level"),
        col_types = "ciciid",
        delim = " "
      ) %>%
        mutate(Sample = str_remove(x, path))
    }
  ))

  index_cpgs <- inner_join(cpg, cpgs)

  cpg_alt <- left_join(index_cpgs, nodes) %>%
    filter(type == "alt") %>%
    arrange(node, pos, strand) %>%
    mutate(base = paste(node, pos, strand, sep = "_")) %>%
    select(base, level, Sample)

  cpg_alt
}

make_cpg_heatmap <- function(cpg_df, plot_path, plot_title) {
  cpg_alt_matrix <- cpg_df %>%
    pivot_wider(id_cols = Sample, names_from = base, values_from = level, values_fill = 0) %>%
    select(-Sample) %>%
    as.matrix() %>%
    t()

    png(plot_path, width = 7, height = 10, units = "in", res = 600)

  lha1 <- rowAnnotation(node = str_split(rownames(cpg_alt_matrix), "_", simplify = T)[, 1], show_legend = F)
  cpg_hm <- Heatmap(cpg_alt_matrix,
    column_labels = NULL, row_labels = NULL, show_column_names = F, show_row_names = F,
    show_row_dend = F, show_column_dend = F, name = "mC level",
    column_title = plot_title,
    row_title = "Methylated nucleotide",
    cluster_rows = F,
    left_annotation = lha1
  )
  draw(cpg_hm)
  dev.off()
  cpg_hm
}

kir <- make_cpg_heatmap(load_alt_cpg("data/bubbles/", "KIR", nodes, cpgs), "plots/kir_alt_mc_full.png", "KIR")
make_cpg_heatmap(load_alt_cpg("data/bubbles/", "HLA-DRB1-DQB1", nodes, cpgs), "plots/hla_drb1_dqb1_alt_mc_full.png", "HLA-DRB1_DQB1")
make_cpg_heatmap(load_alt_cpg("data/bubbles/", "HLA-C-B", nodes, cpgs), "plots/hla_c_b.png", "HLA-C-B")

## sample_cpgs <- lapply(
##   list.files("data/calls/", pattern = "*.nodes", full.names = T),
##   function(x) {
##     snodes <- read_csv(x,
##       col_names = c("node"),
##       )
##     left_join(snodes, cpgs) %>%
##       group_by(node) %>%
##       summarise(ncpg = n()) %>%
##       mutate(sample = str_remove(basename(x), ".fasta.gz.calls.nodes"))
##   }
## )

## sample_cpgs <- bind_rows(sample_cpgs)
## bind_rows(sample_cpgs) %>% write_tsv(sample_cpgs, "data/mg_ncpg_sample.tsv.gz")
sample_cpgs <- read_tsv("data/mg_ncpg_sample.tsv.gz")

ncpgs_sample <- sample_cpgs %>%
  group_by(sample) %>%
  summarise(ncpg_sample = sum(ncpg) / 2)

ncpgs_sample_p <- ggplot(ncpgs_sample) +
  geom_histogram(aes(x = ncpg_sample / 1e6)) +
  labs(
    title = "Number of CpGs in GA4K assemblies",
    x = "Number of CpGs (millions)"
  )
ggsave("plots/ncpgs_per_sample.pdf", ncpgs_sample_p)

ncpgs_sample_alt <- left_join(sample_cpgs, nodes) %>%
  filter(type != "ref") %>%
  group_by(sample) %>%
  summarise(ncpg_sample = sum(ncpg) / 2)

ncpgs_sample_alt_p <- ggplot(ncpgs_sample_alt) +
  geom_histogram(aes(x = ncpg_sample / 1e3)) +
  labs(
    title = "Number of CpGs in GA4K assemblies",
    subtitle = "In non-reference sequences",
    x = "Number of CpGs (thousands)"
  )
ggsave("plots/ncpgs_per_sample_alt.pdf", ncpgs_sample_alt_p)
