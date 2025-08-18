library(ggplot2)
library(tidyverse)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
# outfile_prefix <- args[2]
mapq <- as.numeric(args[2])
cov <- as.numeric(args[3])
annotation <- args[4]

# read annotation
annot <- read.table(annotation, header = TRUE, sep = '\t')
colnames(annot) <- c("accession", "name", "gembasesID", "taxid")
annot$taxid <- as.integer(annot$taxid)
# read paf
a <- read.table(filename, header = TRUE, sep = '\t')
column_names <- c('query_name','query_length','query_start','query_stop','strand','target_name','target_length','target_start','targer_end',
                  'matching_bases','nBases','quality', 'sample','group')
colnames(a) <- column_names

# filter by quality, and extract taxID and accession
a <- a %>% mutate(coverage = (query_stop - query_start) / query_length)
af <- a %>% filter(quality >= mapq) %>% filter(coverage >= cov) %>%
  separate(target_name, into = c("accession", NA, "taxid"),sep = "\\|")
af$taxid <- as.integer(af$taxid)

# join annotation (names)
a.annot <- left_join(af, annot, by = c("accession", "taxid"))
write.table(a.annot, paste0("filtered_hits.tsv"), col.names = TRUE, row.names = FALSE,
    sep = '\t', quote = FALSE)

# number of reads per hit, per sample and group
b <- a.annot %>% group_by(group, sample, taxid, name) %>%
  summarise(counts = n())


# heatmap all
# plot size
n_species <- length(unique(b$name))
n_samples <- length(unique(b$sample))
height_per_species <- 0.4  # inches per species
base_height <- 10           # minimal height in inches
width_per_sample <- 0.7    # inches per sample
base_width <- 12            # minimal width in inches
plot_height <- max(base_height, n_species * height_per_species)
plot_width  <- max(base_width,  n_samples * width_per_sample)

pdf(paste0("all_heatmap.pdf"), width = plot_width, height = plot_height)
ggplot(b %>%  filter(!is.na(name)), aes(x = factor(sample), y = name, fill = log2(counts), label = counts)) +
  geom_tile() +
  geom_text(colour = 'white') +
  facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  theme_classic() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# plot boxplot
pdf(paste0("all_boxplot.pdf"), width = 12, height = plot_height)
ggplot(b %>%  filter(!is.na(name)), aes(x = log2(counts), y = name, colour = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2) +
  theme_bw() +
  labs(y = '')
dev.off()

# heatmap by group
for(gr in setdiff(unique(b$group), "control")){
  print(gr)
  pdf(paste0(gr, "_group_heatmap.pdf"), width = plot_width, height = plot_height)
  p = ggplot(b %>%  filter(!is.na(name)) %>% filter(group %in% c("control", gr)),
         aes(x = factor(sample), y = name, fill = log2(counts), label = counts)) +
    geom_tile() +
    geom_text(colour = 'white') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_bw() +
    labs(x = '', y = '') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(p)
  dev.off()
}

# boxplot by group
for(gr in setdiff(unique(b$group), "control")){
  print(gr)
  pdf(paste0(gr, "_group_boxplot.pdf"), width = plot_width, height = plot_height)
  p = ggplot(b %>%  filter(!is.na(name)) %>% filter(group %in% c("control", gr)),
  aes(x = log2(counts), y = name, colour = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2) +
  theme_bw() +
  labs(y = '')
  print(p)
  dev.off()
}
