set.seed(123)
library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
min_reads <- as.numeric(args[2])
print(min_reads)
# contaminants <- args[2]
# tax_levels <- strsplit(args[3], ",\\s*")[[1]]
# print(tax_levels)

df <- read.table(filename, header = TRUE, sep = '\t')

# summary table of reads by Domain
sum.domains <- df %>%
  mutate(domain = case_when(
    rank == "superkingdom" & taxName == "Viruses" ~ "Viruses",
    rank == "superkingdom" & taxName == "Bacteria" ~ "Bacteria",
    rank == "kingdom" & taxName == "Metazoa" ~ "Human",
    rank == "kingdom"  & taxName == "Fungi" ~ "Fungi",
    rank == "no rank"  & taxName == "unclassified" ~ "Unclassified",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(domain)) %>%
  group_by(sample, domain) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = reads) %>%
  mutate(across(everything(), ~replace_na(., 0)))
write.table(sum.domains, "summary_kingdoms.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# for now, plot at species level
tax_levels <- 'species'
# select rank
for(tax_level in tax_levels) {
  print(paste0("Processing at ", tax_level, " level"))

  # table at TAX_LEVEL, 
  # FILTERS: 
  #        * have at least "min_reads" reads
  #        * exclude human, bradys and others
   a <- df %>% filter(rank == tax_level) %>% select(reads, taxReads, kmers, dup, cov, taxID, taxName, sample, group) %>%
        filter(reads > min_reads) %>%
        filter(!taxName %in% c("Homo sapiens", "synthetic construct"),
        !startsWith(taxName, "Bradyrhizobium"))

  write.table(a, paste0("long_table_at_", tax_level, "_level.tsv"),
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
#   # wide table
  watot <- a %>%
    select(-taxReads, -kmers, -dup, -cov) %>%
    unite(sample_group, sample, group, sep = "_") %>%
    pivot_wider(
      names_from = sample_group,
      values_from = reads,
      values_fill = 0
    )
  write.table(watot, paste0("table_at_", tax_level, "_level_reads.tsv"),
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

#   wadis <- a %>%
#     select(-taxReads) %>%
#     unite(sample_group, sample, group, sep = "_") %>%
#     pivot_wider(
#       names_from = sample_group,
#       values_from = distinctMinimizers,
#       values_fill = 0
#     )
#   write.table(wadis, paste0("table_at_", tax_level, "_level_distinctMinimizers.tsv"),
#         row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

  ## heatmap all species
  # reorder by total reads per species across all samples
  totals <- a %>%
    group_by(taxName) %>%
    summarise(totReads_species = sum(reads))
  a <- a %>% left_join(totals, by = "taxName")
  a$taxName <- factor(a$taxName, levels = totals$taxName[order(totals$totReads_species)])

  # a <- a %>%
  #   mutate(taxName = fct_reorder(taxName, totReads_species))

  # plot size
  n_species <- length(unique(a$taxName))
  n_samples <- length(unique(a$sample))
  height_per_species <- 0.2  # inches per species
  base_height <- 4           # minimal height in inches
  width_per_sample <- 0.5    # inches per sample
  base_width <- 8            # minimal width in inches
  plot_height <- max(base_height, n_species * height_per_species)
  plot_width  <- max(base_width,  n_samples * width_per_sample)

  pdf(paste0(tax_level, "Level_heatmap_reads_fillByCoverage.pdf"), width = plot_width, height = plot_height)
  p = ggplot(a  %>% filter(reads != 0), aes(x = factor(sample), y = taxName, fill = cov, label = reads)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(p)
  dev.off()

#   pdf(paste0(tax_level, "Level_heatmap_distinctMinimizers.pdf"), width = plot_width, height = plot_height)
#   p = ggplot(a %>% filter(distinctMinimizers != 0), aes(x = factor(sample), y = taxName, fill = distinctMinimizers, label = distinctMinimizers)) +
#     geom_tile() +
#     geom_text(colour='white') +
#     labs(x = '', y = '') +
#     facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#   print(p)
#   dev.off()

  # heatmap by group
  for(gr in setdiff(unique(a$group), unique(a$group[grepl("^control", a$group)]))){
    print(gr)
    control_groups <- unique(a$group[grepl("^control", a$group)])
    pdf(paste0(gr, "_", tax_level, "Level_heatmap_reads_fillByCoverage.pdf"), height = plot_height, width = plot_width)
    p = ggplot(a %>% filter(reads != 0) %>% filter(group %in% c(gr, control_groups)), 
          aes(x = factor(sample), y = taxName, fill = cov, label = reads)) +
      geom_tile() +
      geom_text(colour='white') +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()
  }

  # reads vs kmers for each sample
  for(s in unique(a$sample)){
    print(s)
    pdf(paste0(s, "_", tax_level, "Level_reads_vs_kmers.pdf"), 
    width = 10, height = 10)
    p = ggplot(a %>% filter(sample == s), 
          aes(x = log2(reads), y = log2(kmers), label = taxName, size = cov)) +
      geom_point() +
      geom_text_repel() +
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(limits = c(0, NA)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 14)) +
      ggtitle(s)
    print(p)
    dev.off()
  }

#   ### remove contaminants ###
#   if (contaminants != "") {
#     print('Contaminant file provided')
#     conta <- read.table(contaminants, header = FALSE, sep = '\t')$V1
#     b <- a %>% filter(!taxid %in% conta)

#     # plot size
#     n_species <- length(unique(a$taxName))
#     n_samples <- length(unique(a$sample))
#     height_per_species <- 0.2  # inches per species
#     base_height <- 4           # minimal height in inches
#     width_per_sample <- 0.5    # inches per sample
#     base_width <- 8            # minimal width in inches
#     plot_height <- max(base_height, n_species * height_per_species)
#     plot_width  <- max(base_width,  n_samples * width_per_sample)

#     pdf(paste0(tax_level, "Level_heatmap_taxReads_contaminantsRemoved.pdf"), width = plot_width, height = plot_height)
#     ggplot(b  %>% filter(taxReads != 0), aes(x = factor(sample), y = taxName, fill = taxReads, label = taxReads)) +
#       geom_tile() +
#       geom_text(colour='white') +
#       labs(x = '', y = '') +
#       facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#     dev.off()

#     pdf(paste0(tax_level, "Level_heatmap_distinctMinimizers_contaminantsRemoved.pdf"), width = plot_width, height = plot_height)
#     ggplot(b %>% filter(distinctMinimizers != 0), aes(x = factor(sample), y = taxName, fill = distinctMinimizers, label = distinctMinimizers)) +
#       geom_tile() +
#       geom_text(colour='white') +
#       labs(x = '', y = '') +
#       facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#     dev.off()

#     # heatmap by group
#     for(gr in setdiff(unique(a$group), unique(a$group[grepl("^control", a$group)]))){
#       print(gr)
#       control_groups <- unique(a$group[grepl("^control", a$group)])
#       pdf(paste0(gr, "_", tax_level, "Level_heatmap_taxReads_contaminantsRemoved.pdf"), height = plot_height, width = plot_width)
#       p = ggplot(b %>% filter(taxReads != 0) %>% filter(group %in% c(gr, control_groups)), 
#             aes(x = factor(sample), y = taxName, fill = taxReads, label = taxReads)) +
#         geom_tile() +
#         geom_text(colour='white') +
#         labs(x = '', y = '') +
#         facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
#         theme_classic() +
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#       print(p)
#       dev.off()
#     }

#     # total counts vs distinct minimizers for each sample
#     for(s in unique(a$sample)){
#       print(s)
#       pdf(paste0(s, "_", tax_level, "Level_taxReads_vs_distinctMinimizers_contaminantsRemoved.pdf"), 
#       width = 10, height = 10)
#       p = ggplot(b %>% filter(sample == s), 
#             aes(x = log2(taxReads+1), y = log2(distinctMinimizers+1), label = taxName)) +
#         geom_point() +
#         geom_text_repel() +
#         scale_y_continuous(limits = c(0, NA)) +
#         scale_x_continuous(limits = c(0, NA)) +
#         theme_classic() +
#         theme(axis.title = element_text(size = 14),
#               axis.text = element_text(size = 14)) +
#         ggtitle(s)
#       print(p)
#       dev.off()
#     }
#   } else {
#     print('No contaminant file provided')
#   }
 }
