#!/usr/bin/env Rscript
set.seed(123)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------------
# Parse command-line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: script.R <input.tsv>")
}
meta <- args[1]
level <- args[2]
file_list <- args[3:length(args)]


# ---------------------------
#  METADATA
# ---------------------------
metadata = read.table(meta, header = TRUE, sep =',')
metadata = metadata %>% select(sample, group)
metadata$sample = as.character(metadata$sample)

# ---------------------------
# COMBINED TABLES
# ---------------------------
df = data.frame()

for(f in file_list){
    s = sub("_.*", "", f)
    a = read.table(f, header = TRUE, sep = '\t')
    a = a %>% mutate(sample = s)
    df = rbind(df, a)
}
df$sample = as.character(df$sample)
dfg = left_join(df, metadata, by = "sample")

write.table(dfg, paste0('combined_table_by_', level, '.tsv'),
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)

# ---------------------------
# PLOTS FROM TAXID TABLE
# ---------------------------
if(level == 'taxid') {
    ## by species
    # reorder by total
    totals <- dfg %>%
        group_by(species) %>%
        summarise(totReads_species = sum(nBases_covered))
    dfg <- dfg %>% left_join(totals, by = "species")
    dfg$species <- factor(dfg$species, levels = totals$species[order(totals$totReads_species)])

    # plot size
    n_species <- length(unique(dfg$species))
    n_samples <- length(unique(dfg$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    # HETMAP SPECIES BY SAMPLE/GROUP: BASES COVERED
    pdf(paste0('heatmap_nBasesCovered_by_', level, '.pdf'), height = plot_height, width = plot_width)
    p = ggplot(dfg %>% filter(!is.na(species)), aes(x = sample, y = species, fill = nBases_covered, label = nBases_covered)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()


    ## by genus
    dfgg <- dfg %>% group_by(sample, genus) %>%
    summarise(nBases_covered = sum(nBases_covered),
                group = first(group))

    # reorder by total
    totalsg <- dfgg %>%
        group_by(genus) %>%
        summarise(totReads_genus = sum(nBases_covered))
    dfgg <- dfgg %>% left_join(totalsg, by = "genus")
    dfgg$genus <- factor(dfgg$genus, levels = totalsg$genus[order(totalsg$totReads_genus)])

    # plot size
    n_species <- length(unique(dfg$genus))
    n_samples <- length(unique(dfg$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    # HETMAP GENUS BY SAMPLE/GROUP: BASES COVERED
    pdf(paste0('heatmap_nBasesCovered_by_', level, '_genusLevel.pdf'), height = plot_height, width = plot_width)
    p = ggplot(dfgg %>% filter(!is.na(genus)), aes(x = sample, y = genus, fill = nBases_covered, label = nBases_covered)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()

}