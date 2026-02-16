set.seed(123)
library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
min_reads <- as.numeric(args[2])
print(paste0("Taxa with less than ", min_reads, "were removed from plots"))

df <- read.table(filename, header = TRUE, sep = '\t', quote = "", comment.char="")

# wide table of reads per sample and species
w.df <- df %>% filter(rank == 'species') %>% 
  select(reads, sample, group, taxID, taxName) %>%
  unite(sample_group, sample, group, sep = "_") %>%
  pivot_wider(
    names_from = sample_group,
    values_from = reads,
    values_fill = 0
  )
write.table(w.df, paste0("krakenuniq_read_table_at_species_level.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# summary table of reads by Domain
sum.domains <- df %>%
  mutate(domain = case_when(
    rank == "superkingdom" & taxName == "Viruses" ~ "krakenuniq_Viruses",
    rank == "superkingdom" & taxName == "Bacteria" ~ "krakenuniq_Bacteria",
    rank == "kingdom" & taxName == "Metazoa" ~ "krakenuniq_Human",
    rank == "kingdom"  & taxName == "Fungi" ~ "krakenuniq_Fungi",
    rank == "no rank"  & taxName == "unclassified" ~ "krakenuniq_Unclassified",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(domain)) %>%
  group_by(sample, domain) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = reads) %>%
  mutate(across(everything(), ~replace_na(., 0)))
write.table(sum.domains, "krakenuniq_summary_kingdoms.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

# Remove human from every plot and table
a.hsRem <- df %>%
    filter(taxName !=  "Homo sapiens")

# FILTERS: 
#        * have at least "min_reads" reads
#        * exclude bradys (and others ?) manual removal here; a nextflow.config param contaminants still not implemented
bef.reads <- a.hsRem %>% 
    filter(rank == 'species') %>%
    group_by(sample) %>%
    summarise(beforeReads = sum(reads),
              krakenuniq_beforeTaxa = n() + 1) # add 1 for human

a <- a.hsRem %>% filter(rank == 'species') %>% 
    select(reads, taxReads, kmers, dup, cov, taxID, taxName, sample, group) %>%
    filter(reads > min_reads) %>%
    filter(!taxName %in% c("synthetic construct"),
    !startsWith(taxName, "Bradyrhizobium"))
# write.table(a, paste0("long_table_at_species_level.tsv"),
#         row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

aft.reads <- a %>%
  group_by(sample) %>%
  summarise(afterReads = sum(reads),
            krakenuniq_afterTaxa = n())

rem.reads <- left_join(bef.reads, aft.reads, by = "sample") %>%
    mutate(krakenuniq_removedReadsFromPlots = beforeReads - afterReads)
write.table(rem.reads, paste0("krakenuniq_removedReadsFromPlots.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

# wide table of reads per sample and species
watot <- a %>%
  select(-taxReads, -kmers, -dup, -cov) %>%
  unite(sample_group, sample, group, sep = "_") %>%
  pivot_wider(
    names_from = sample_group,
    values_from = reads,
    values_fill = 0
  )
write.table(watot, paste0("krakenuniq_read_table_at_species_level_postCleaning.tsv"),
      row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

## heatmap all species
# reorder by total reads per species across all samples
totals <- a %>%
  group_by(taxName) %>%
  summarise(totReads_species = sum(reads))
a <- a %>% left_join(totals, by = "taxName")
a$taxName <- factor(a$taxName, levels = totals$taxName[order(totals$totReads_species)])

# plot size
n_species <- length(unique(a$taxName))
n_samples <- length(unique(a$sample))
height_per_species <- 0.2  # inches per species
base_height <- 4           # minimal height in inches
width_per_sample <- 0.5    # inches per sample
base_width <- 8            # minimal width in inches
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

pdf(paste0("heatmap_reads_fillByCoverage.pdf"), width = plot_width, height = plot_height)
p = ggplot(a  %>% filter(reads != 0), aes(x = factor(sample), y = taxName, fill = cov, label = reads)) +
  geom_tile() +
  geom_text(colour='white') +
  labs(x = '', y = '') +
  facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p)
dev.off()

# heatmap by group
for(gr in setdiff(unique(a$group), unique(a$group[grepl("^control", a$group)]))){
  print(gr)
  control_groups <- unique(a$group[grepl("^control", a$group)])
  pdf(paste0(gr, "_heatmap_reads_fillByCoverage.pdf"), height = plot_height, width = plot_width)
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
  pdf(paste0(s, "_reads_vs_kmers.pdf"), 
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
