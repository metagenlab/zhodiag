library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outfile_prefix <- args[2]

a <- read.table(filename, header = TRUE, sep = '\t')

# select rank species
a <- a %>% filter(rank == "S")

# heatmap all
# plot size
n_species <- length(unique(a$taxonomy))
n_samples <- length(unique(a$sample))
height_per_species <- 0.2  # inches per species
base_height <- 4           # minimal height in inches
width_per_sample <- 0.5    # inches per sample
base_width <- 8            # minimal width in inches
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

pdf(paste0(outfile_prefix, "_heatmap_all.pdf"), width = plot_width, height = plot_height)
ggplot(a, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

### remove homo sapiens ###
b <- a %>% filter(taxonomy != 'Homo sapiens')

# heatmap without human
pdf(paste0(outfile_prefix, "_heatmap_exclHuman.pdf"), width = plot_width, height = plot_height)
ggplot(b, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# by group facet
pdf(paste0(outfile_prefix, "_heatmap_exclHuman_facetGroups.pdf"), width = plot_width, height = plot_height)
ggplot(b, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# boxplot totalCounts
b$group <- factor(b$group)
pdf(paste0(outfile_prefix, "_boxplot_totalCounts_exclHuman.pdf"), width = 12, height = plot_height)
ggplot(b, aes(x = taxonomy, y = totalCounts, colour = group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "", y = "Total Counts", colour = "Group") +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") + 
  # guides(colour = "none") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
dev.off()

# boxplot distinctMinimizers
pdf(paste0(outfile_prefix, "_boxplot_distinctMinimizers_exclHuman.pdf"), width = 12, height = plot_height)
ggplot(b, aes(x = taxonomy, y = distinctMinimizers, colour = group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "", y = "Distinct Minimizers", colour = "Group") +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") + 
  # guides(colour = "none") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
dev.off()

# total counts vs distinct minimizers for each sample
for(s in unique(b$sample)){
#   print(s)
    pdf(paste0(outfile_prefix, "_", s, "_totalCounts_vs_distinctMinimizers.pdf"), 
    width = 10, height = 10)
  p = ggplot(b %>% filter(sample == s), 
         aes(x = totalCounts, y = distinctMinimizers, label = taxonomy)) +
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

#### FILTER ###
# remove what is found only in controls
taxonomy_summary <- b %>%
  group_by(taxonomy, sample) %>%
  summarise(sample_sum = sum(totalCounts), .groups = "drop")
taxonomy_wide <- taxonomy_summary %>%
  tidyr::pivot_wider(names_from = sample, values_from = sample_sum, values_fill = 0)
control_cols <- b %>% filter(group == 'control') %>% pull(sample) %>% unique()
non_control_cols <- setdiff(names(taxonomy_wide), c(control_cols, "taxonomy"))
taxa_to_remove <- taxonomy_wide %>%
  filter(control_cols > 0, rowSums(across(all_of(non_control_cols))) == 0) %>%
  pull(taxonomy)

bf <- b %>%
  filter(!taxonomy %in% taxa_to_remove)

n_species <- length(unique(bf$taxonomy))
n_samples <- length(unique(bf$sample))
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

# heatmap together
pdf(paste0(outfile_prefix, "_heatmap_exclControls.pdf"), width = plot_width, height = plot_height)
ggplot(bf, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# heatmap facet
pdf(paste0(outfile_prefix, "_heatmap_exclControls_facetGroups.pdf"), width = plot_width, height = plot_height)
ggplot(bf, aes(x = factor(sample), y = taxonomy, fill = totalCounts, label = totalCounts)) +
  geom_tile() +
  geom_text(colour='white') +
  facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# # boxplot totalCounts
# bf$group <- factor(b$group)
# pdf(paste0(outfile_prefix, "_boxplot_totalCounts_exclControls.pdf"), width = 12, height = plot_height)
# ggplot(bf, aes(x = taxonomy, y = totalCounts, fill = group, colour = group)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +
#   labs(x = "", y = "Total Counts", fill = "Group") +
#   theme_bw() +
#   guides(colour = "none") +
#   # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   coord_flip()
# dev.off()

# # boxplot distinctMinimizers
# pdf(paste0(outfile_prefix, "_boxplot_distinctMinimizers_exclControls.pdf"), width = 12, height = plot_height)
# ggplot(bf, aes(x = taxonomy, y = distinctMinimizers, fill = group, colour = group)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +
#   labs(x = "", y = "Distinct Minimizers", fill = "Group") +
#   theme_bw() +
#   guides(colour = "none") +
#   # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   coord_flip()
# dev.off()

# total counts vs distinct minimizers for each sample
for(s in unique(b$sample)){
    pdf(paste0(outfile_prefix, "_", s, "_totalCounts_vs_distinctMinimizers_exclControls.pdf"), 
    width = 10, height = 10)
  p = ggplot(b %>% filter(sample == s), 
         aes(x = totalCounts, y = distinctMinimizers, label = taxonomy)) +
    geom_point() +
    geom_text_repel() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ggtitle(s)
  print(p)
}
