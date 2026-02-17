set.seed(123)
library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
min_reads <- as.numeric(args[2])
contaminants <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL

print(paste0("Taxa with less than ", min_reads, " reads were removed from plots"))

df <- read.table(filename, header = TRUE, sep = '\t', quote = "", comment.char="")

if (!is.null(contaminants)) {
  contams <- read.table(contaminants, header = FALSE, sep = '\t')
  colnames(contams) <- c('contaminant', 'code')
  print("The following taxa were removed: ")
  contams
  write.table(contams, "removed_contaminants.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')
}

# wide table of reads per sample and species
w.df <- df %>% filter(rank == 'S') %>% 
  select(totalCounts, sample, group, taxid, taxonomy) %>%
  unite(sample_group, sample, group, sep = "_") %>%
  pivot_wider(
    names_from = sample_group,
    values_from = totalCounts,
    values_fill = 0
  )
write.table(w.df, paste0("kraken2_read_table_at_species_level.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

# summary table of nReads by Domain
sum.domains <- df %>%
  mutate(domain = case_when(
    rank == "R1" & taxonomy == "Viruses" ~ "kraken2_Viruses",
    rank == "R2" & taxonomy == "Bacteria" ~ "kraken2_Bacteria",
    rank == "K" & taxonomy == "Metazoa" ~ "kraken2_Human",
    rank == "K"  & taxonomy == "Fungi" ~ "kraken2_Fungi",
    rank == "U"  & taxonomy == "unclassified" ~ "kraken2_Unclassified",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(domain)) %>%
  group_by(sample, domain) %>%
  summarise(nreads = sum(totalCounts), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = nreads) %>%
  mutate(across(everything(), ~replace_na(., 0)))
write.table(sum.domains, "kraken2_summary_kingdoms.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# Remove human from every plot and table
a.hsRem <- df %>%
    filter(taxonomy !=  "Homo sapiens")

# FILTERS: 
#        * have at least "min_reads" reads
#        * exclude taxa specified in params.contaminants file
#        * only detected in control group
bef.reads <- a.hsRem %>% 
    filter(rank == 'S') %>%
    group_by(sample) %>%
    summarise(beforeReads = sum(totalCounts),
              kraken2_beforeTaxa = n() + 1) # add 1 for human

a.filt1 <- a.hsRem %>% filter(rank == 'S') %>% 
    select(totalCounts, distinctMinimizers, taxid, taxonomy, sample, group) %>%
    filter(totalCounts > min_reads)

if (!is.null(contaminants)) {
  sp.filter = contams %>% filter(code == 'species') %>% pull(contaminant)
  gn.filter = contams %>% filter(code == 'genus') %>% pull(contaminant)
  tax.filter = contams %>% filter(code == 'taxid') %>% pull(contaminant)

  a.filt1 <- a.filt1 %>% 
      filter(!taxonomy %in% sp.filter,
              !Reduce(`|`, lapply(gn.filter, function(g) startsWith(taxonomy, g))),
              !taxid %in% tax.filter)
}

a <- a.filt1 %>%
  group_by(taxonomy) %>%
  mutate(
    control_reads = sum(totalCounts[group == "control"]),
    other_reads   = sum(totalCounts[group != "control"])
  ) %>%
  ungroup() %>%
  filter(!(control_reads > 0 & other_reads == 0)) %>%
  select(-control_reads, -other_reads)
# write.table(a, paste0("long_table_at_species_level.tsv"),
#         row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

aft.reads <- a %>%
  group_by(sample) %>%
  summarise(afterReads = sum(totalCounts),
            kraken2_afterTaxa = n())

rem.reads <- left_join(bef.reads, aft.reads, by = "sample") %>%
    mutate(kraken2_removedReadsFromPlots = beforeReads - afterReads)
write.table(rem.reads, paste0("kraken2_removedReadsFromPlots.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

# wide table of reads per sample and species
watot <- a %>%
  select(-distinctMinimizers) %>%
  unite(sample_group, sample, group, sep = "_") %>%
  pivot_wider(
    names_from = sample_group,
    values_from = totalCounts,
    values_fill = 0
  )
write.table(watot, paste0("kraken2_read_table_at_species_level_postCleaning.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


## heatmap all species
# reorder by total reads per species across all samples
totals <- a %>%
  group_by(taxonomy) %>%
  summarise(totReads_species = sum(totalCounts))
a <- a %>% left_join(totals, by = "taxonomy")
a$taxonomy <- factor(a$taxonomy, levels = totals$taxonomy[order(totals$totReads_species)])

# plot size
height_per_species <- 0.2  # inches per species
base_height <- 4           # minimal height in inches
width_per_sample <- 0.5    # inches per sample
base_width <- 8            # minimal width in inches

n_species <- length(unique(a$taxonomy))
n_samples <- length(unique(a$sample))
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

pdf("heatmap_totalCounts.pdf", width = plot_width, height = plot_height)
p = ggplot(a, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  labs(x = '', y = '') +
  facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p)
dev.off()

#   pdf(paste0(tax_level, "Level_heatmap_distinctMinimizers.pdf"), width = plot_width, height = plot_height)
#   p = ggplot(a %>% filter(distinctMinimizers != 0), aes(x = factor(sample), y = taxonomy, fill = distinctMinimizers, label = distinctMinimizers)) +
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
  dtpm <- a %>% filter(group %in% c(gr, control_groups))

  # filter out taxa only in control vs group
  dtp <- dtpm %>%
  group_by(taxonomy) %>%
  mutate(
    control_reads = sum(totalCounts[group == "control"]),
    other_reads   = sum(totalCounts[group != "control"])
  ) %>%
  ungroup() %>%
  filter(!(control_reads > 0 & other_reads == 0)) %>%
  select(-control_reads, -other_reads)

  # reorder by total reads per species across all samples
  totals.gr <- dtp %>%
    group_by(taxonomy) %>%
    summarise(totReads_species = sum(totalCounts))
  dtp <- dtp %>% left_join(totals.gr, by = "taxonomy")
  dtp$taxonomy <- factor(dtp$taxonomy, levels = totals.gr$taxonomy[order(totals.gr$totReads_species)])

  # plot size
  n_species <- length(unique(dtp$taxonomy))
  n_samples <- length(unique(dtp$sample))
  plot_height <- max(base_height, n_species * height_per_species)
  plot_width  <- max(base_width,  n_samples * width_per_sample)

  pdf(paste0(gr, "_heatmap_totalCounts.pdf"), height = plot_height, width = plot_width)
  p = ggplot(dtp, 
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
  pdf(paste0(s, "_totalCounts_vs_distinctMinimizers.pdf"), 
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
