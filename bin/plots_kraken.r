library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
# outfile_prefix <- args[2]
# contaminants <- as.numeric(strsplit(args[2], ",\\s*")[[1]])
contaminants <- args[2]
# tax_level <- as.character(args[3])
tax_levels <- strsplit(args[3], ",\\s*")[[1]]
print(tax_levels)

df <- read.table(filename, header = TRUE, sep = '\t')

# summary table of nReads by Domain
sum.domains <- df %>%
  mutate(domain = case_when(
    rank == "R1" & taxonomy == "Viruses" ~ "Viruses",
    rank == "R2" & taxonomy == "Bacteria" ~ "Bacteria",
    rank == "K" & taxonomy == "Metazoa" ~ "Human",
    rank == "K"  & taxonomy == "Fungi" ~ "Fungi",
    rank == "U"  & taxonomy == "unclassified" ~ "Unclassified",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(domain)) %>%
  group_by(sample, domain) %>%
  summarise(nreads = sum(totalCounts), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = nreads) %>%
  mutate(across(everything(), ~replace_na(., 0)))
write.table(sum.domains, "summary_kingdoms.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# select rank
for(tax_level in tax_levels) {
  print(paste0("Processing at ", tax_level, " level"))
  selected_rank <- toupper(substr(tax_level, 1, 1))
  a <- df %>% filter(rank == selected_rank) %>% select(totalCounts, distinctMinimizers, taxid, taxonomy, sample, group)

  write.table(a, paste0("long_table_at_", tax_level, "_level.tsv"),
        row.names = FALSE, col.names = TRUE, sep = '\t', ,quote = FALSE)
  # wide table
  watot <- a %>%
    select(-distinctMinimizers) %>%
    unite(sample_group, sample, group, sep = "_") %>%
    pivot_wider(
      names_from = sample_group,
      values_from = totalCounts,
      values_fill = 0
    )
  write.table(watot, paste0("table_at_", tax_level, "_level_totalCounts.tsv"),
        row.names = FALSE, col.names = TRUE, sep = '\t', ,quote = FALSE)

  wadis <- a %>%
    select(-totalCounts) %>%
    unite(sample_group, sample, group, sep = "_") %>%
    pivot_wider(
      names_from = sample_group,
      values_from = distinctMinimizers,
      values_fill = 0
    )
  write.table(wadis, paste0("table_at_", tax_level, "_level_distinctMinimizers.tsv"),
        row.names = FALSE, col.names = TRUE, sep = '\t', ,quote = FALSE)

  ## heatmap all species
  # reorder by total reads per species across all samples
  totals <- a %>%
    group_by(taxonomy) %>%
    summarise(totReads_species = sum(totalCounts))
  a <- a %>% left_join(totals, by = "taxonomy")
  a$taxonomy <- factor(a$taxonomy, levels = totals$taxonomy[order(totals$totReads_species)])

  # a <- a %>%
  #   mutate(taxonomy = fct_reorder(taxonomy, totReads_species))

  # plot size
  n_species <- length(unique(a$taxonomy))
  n_samples <- length(unique(a$sample))
  height_per_species <- 0.2  # inches per species
  base_height <- 4           # minimal height in inches
  width_per_sample <- 0.5    # inches per sample
  base_width <- 8            # minimal width in inches
  plot_height <- max(base_height, n_species * height_per_species)
  plot_width  <- max(base_width,  n_samples * width_per_sample)

  pdf(paste0(tax_level, "Level_heatmap_totalCounts.pdf"), width = plot_width, height = plot_height)
  ggplot(a  %>% filter(totalCounts != 0), aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  dev.off()

  pdf(paste0(tax_level, "Level_heatmap_distinctMinimizers.pdf"), width = plot_width, height = plot_height)
  ggplot(a %>% filter(distinctMinimizers != 0), aes(x = factor(sample), y = taxonomy, fill = distinctMinimizers, label = distinctMinimizers)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  dev.off()

  # heatmap by group
  for(gr in setdiff(unique(a$group), 'control')){
    print(gr)
    pdf(paste0(gr, "_", tax_level, "Level_heatmap_totalCounts.pdf"), height = plot_height, width = plot_width)
    p = ggplot(a %>% filter(totalCounts != 0) %>% filter(group %in% c(gr, "control")), 
          aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
      geom_tile() +
      geom_text(colour='white') +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()
  }

  # total counts vs distinct minimizers for each sample
  for(s in unique(a$sample)){
    print(s)
    pdf(paste0(s, "_", tax_level, "Level_totalCounts_vs_distinctMinimizers.pdf"), 
    width = 10, height = 10)
    p = ggplot(a %>% filter(sample == s), 
          aes(x = log2(totalCounts+1), y = log2(distinctMinimizers+1), label = taxonomy)) +
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

  ### remove contaminants ###
  if (contaminants != "") {
    print('Contaminant file provided')
    conta <- read.table(contaminants, header = FALSE, sep = '\t')$V1
    b <- a %>% filter(!taxid %in% conta)

    # plot size
    n_species <- length(unique(a$taxonomy))
    n_samples <- length(unique(a$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    pdf(paste0(tax_level, "Level_heatmap_totalCounts_contaminantsRemoved.pdf"), width = plot_width, height = plot_height)
    ggplot(b  %>% filter(totalCounts != 0), aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
      geom_tile() +
      geom_text(colour='white') +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    dev.off()

    pdf(paste0(tax_level, "Level_heatmap_distinctMinimizers_contaminantsRemoved.pdf"), width = plot_width, height = plot_height)
    ggplot(b %>% filter(distinctMinimizers != 0), aes(x = factor(sample), y = taxonomy, fill = distinctMinimizers, label = distinctMinimizers)) +
      geom_tile() +
      geom_text(colour='white') +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    dev.off()

    # heatmap by group
    for(gr in setdiff(unique(a$group), 'control')){
      print(gr)
      pdf(paste0(gr, "_", tax_level, "Level_heatmap_totalCounts_contaminantsRemoved.pdf"), height = plot_height, width = plot_width)
      p = ggplot(b %>% filter(totalCounts != 0) %>% filter(group %in% c(gr, "control")), 
            aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
        geom_tile() +
        geom_text(colour='white') +
        labs(x = '', y = '') +
        facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      print(p)
      dev.off()
    }

    # total counts vs distinct minimizers for each sample
    for(s in unique(a$sample)){
      print(s)
      pdf(paste0(s, "_", tax_level, "Level_totalCounts_vs_distinctMinimizers_contaminantsRemoved.pdf"), 
      width = 10, height = 10)
      p = ggplot(b %>% filter(sample == s), 
            aes(x = log2(totalCounts+1), y = log2(distinctMinimizers+1), label = taxonomy)) +
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
  } else {
    print('No contaminant file provided')
  }
}
