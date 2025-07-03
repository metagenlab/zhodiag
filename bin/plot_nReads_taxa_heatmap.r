library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outfile_prefix <- args[2]

f <- read.table(filename, header = TRUE, sep = '\t')
tax_cols <- c("taxid", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
f <- f %>%
  filter(if_all(all_of(tax_cols), ~ !is.na(.) & . != ""))
# add total reads; append taxID to make unique rows, sort by tot reads; get top 100
f <- f %>% mutate(totalReads = rowSums(f[,10:ncol(f)]))
f <- f %>% mutate(species_label = paste0(species, " (", taxid, ")"))
f <- f %>% arrange(desc(totalReads))
f$species_label <- factor(f$species_label, levels = rev(f$species_label))
if (nrow(f) > 100){
  f <- f[1:100,]
}

f$species <- f$species_label
fl <- pivot_longer(f, cols = 10:(ncol(f) - 1), names_to = 'sample', values_to = 'nReads')

# normalize by sample
fl <- fl %>%
  group_by(sample) %>%
  mutate(normReads = nReads / max(nReads, na.rm = TRUE)) %>%
  ungroup()


# plot size
n_species <- length(unique(fl$species))
n_samples <- length(unique(fl$sample))
height_per_species <- 0.2  # inches per species
base_height <- 4           # minimal height in inches
width_per_sample <- 0.7    # inches per sample
base_width <- 8            # minimal width in inches
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

# all taxa
pdf(paste0(outfile_prefix, "_heatmap_nTop100.pdf"), width = plot_width, height = plot_height)
ggplot(fl, aes(x = sample, y = species, fill = normReads, label = nReads)) +
  geom_tile() +
  geom_text(size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue", guide = FALSE) +
  theme_minimal() +
  ggtitle("Top hits") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(fill = "n Reads")
dev.off()

#exclude human, unclassified, ambiguous, bradyrhizobium
pdf(paste0(outfile_prefix, "_heatmap_nTop100_exclHumanAmbigUnclasBrady.pdf"), width = plot_width, height = plot_height)
fl_ex <- fl %>% 
  filter(!str_detect(species, "^(Bradyrhizobium|Homo|ambiguous|unclassified)"))

fl_ex <- fl_ex %>%
  group_by(sample) %>%
  mutate(normReads = nReads / max(nReads, na.rm = TRUE)) %>%
  ungroup()
  
ggplot(fl_ex, aes(x = sample, y = species, fill = normReads, label = nReads)) +
  geom_tile() +
  geom_text(size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue", guide = FALSE) +
  theme_minimal() +
  ggtitle("Top hits", subtitle = "excluding human, ambiguous, unclassified and Bradyrhizobiums") +
  labs(fill = "n Reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()