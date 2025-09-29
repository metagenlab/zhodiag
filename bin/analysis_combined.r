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
# COMBINE TABLES
# ---------------------------
df = data.frame()

for(f in file_list){
    s = sub("_summary.*", "", f)
    print(s)
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
cov_colors <- c(
  "0–300"      = "#f7fbff",
  "300–500"    = "#deebf7",
  "500–1000"   = "#c6dbef",
  "1000–2000"  = "#9ecae1",
  "2000-5000"  = "#6baed6",
  "5000-10000" = "#3182bd",
  "10000+"     = "#08519c"
)

if(level == 'taxid') {
    ## by species
    # reorder by total

  dfg.hm <- dfg %>%
    group_by(sample, species) %>%
    mutate(mappedReads = sum(mappedReads, na.rm = TRUE),
          nBases_covered = sum(nBases_covered)) %>%
    ungroup() %>%
    filter(!is.na(species)) %>%
    filter(mappedReads > 1, nBases_covered > 151) %>%
    mutate(species = reorder(species, mappedReads)) %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
      labels = c("0–300", "300–500", "500–1000", "1000–2000", "2000-5000", "5000-10000", "10000+"),
      right = FALSE
    ))

    # totals <- dfg %>%
    #     group_by(species) %>%
    #     summarise(totReads_species = sum(nBases_covered))
    # dfg <- dfg %>% left_join(totals, by = "species")
    # dfg$species <- factor(dfg$species, levels = totals$species[order(totals$totReads_species)])

    # plot size
    n_species <- length(unique(dfg.hm$species))
    n_samples <- length(unique(dfg.hm$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 10           # minimal height in inches
    width_per_sample <- 1    # inches per sample
    base_width <- 10            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    # HETMAP SPECIES: BASES COVERED and MAPPED READS
    pdf(paste0('heatmap_by_', level, '_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
    p = ggplot(dfg.hm, 
    aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour='white') +
    scale_fill_manual(values = cov_colors) +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()


    ## by genus
    dfgg.hm <- dfg %>%
    group_by(sample, genus) %>%
    mutate(mappedReads = sum(mappedReads, na.rm = TRUE),
            nBases_covered = sum(nBases_covered)) %>%
    filter(!is.na(genus)) %>%
    filter(mappedReads > 1, nBases_covered > 151) %>%
    ungroup() %>%
    mutate(genus = reorder(genus, mappedReads)) %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
      labels = c("0–300", "300–500", "500–1000", "1000–2000", "2000-5000", "5000-10000", "10000+"),
      right = FALSE
    ))


    # dfgg <- dfg %>% group_by(sample, genus) %>%
    # summarise(nBases_covered = sum(nBases_covered),
    #             group = first(group))

    # # reorder by total
    # totalsg <- dfgg %>%
    #     group_by(genus) %>%
    #     summarise(totReads_genus = sum(nBases_covered))
    # dfgg <- dfgg %>% left_join(totalsg, by = "genus")
    # dfgg$genus <- factor(dfgg$genus, levels = totalsg$genus[order(totalsg$totReads_genus)])

    # plot size
    n_species <- length(unique(dfgg.hm$genus))
    n_samples <- length(unique(dfgg.hm$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    # HETMAP GENUS BY SAMPLE/GROUP: BASES COVERED
    pdf(paste0('heatmap_nBasesCovered_by_', level, '_genusLevel.pdf'), height = plot_height, width = plot_width)
    p = ggplot(dfgg.hm %>% filter(!is.na(genus)), aes(x = sample, y = genus, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour='white') +
    labs(x = '', y = '') +
    scale_fill_manual(values = cov_colors) +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()


  # VIRUS
    virus <- dfg %>% filter(domain == "Viruses") %>%
    filter(!is.na(species)) %>%
    filter(mappedReads > 1, nBases_covered > 151) %>%
    mutate(species = reorder(species, mappedReads)) %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
      labels = c("0–300", "300–500", "500–1000", "1000–2000", "2000-5000", "5000-10000", "10000+"),
      right = FALSE
    ))



    # plot size
    n_species <- length(unique(virus$species))
    n_samples <- length(unique(virus$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)


    pdf(paste0('heatmap_VIRUSES_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
    v.plot = ggplot(virus , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour = 'white') +
    facet_grid(~group, scales = 'free_x', space = 'free') +
    labs(x = '', y = '') +
    theme_classic() +
    scale_fill_manual(values = cov_colors) +
    theme(axis.text.y = element_text(size =12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(v.plot)
    dev.off()


  # EUKARYOTA
    euka <- dfg %>% filter(domain == "Eukaryota") %>%
    filter(!is.na(species)) %>%
    filter(mappedReads > 1, nBases_covered > 151) %>%
    mutate(species = reorder(species, mappedReads)) %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
      labels = c("0–300", "300–500", "500–1000", "1000–2000", "2000-5000", "5000-10000", "10000+"),
      right = FALSE
    ))



    # plot size
    n_species <- length(unique(euka$species))
    n_samples <- length(unique(euka$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)


    pdf(paste0('heatmap_EUKARYOTA_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
    e.plot = ggplot(euka , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour = 'white') +
    facet_grid(~group, scales = 'free_x', space = 'free') +
    labs(x = '', y = '') +
    scale_fill_manual(values = cov_colors) +
    theme_classic() +
    theme(axis.text.y = element_text(size =12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(e.plot)
    dev.off()


  # EUKARYOTA
    bact <- dfg %>% filter(domain == "Bacteria") %>%
    filter(!is.na(species)) %>%
    filter(mappedReads > 1, nBases_covered > 151) %>%
    mutate(species = reorder(species, mappedReads)) %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
      labels = c("0–300", "300–500", "500–1000", "1000–2000", "2000-5000", "5000-10000", "10000+"),
      right = FALSE
    ))



    # plot size
    n_species <- length(unique(bact$species))
    n_samples <- length(unique(bact$sample))
    height_per_species <- 0.2  # inches per species
    base_height <- 4           # minimal height in inches
    width_per_sample <- 0.5    # inches per sample
    base_width <- 8            # minimal width in inches
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)


    pdf(paste0('heatmap_BACTERIA_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
    b.plot = ggplot(bact , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour = 'white') +
    facet_grid(~group, scales = 'free_x', space = 'free') +
    labs(x = '', y = '') +
    scale_fill_manual(values = cov_colors) +
    theme_classic() +
    theme(axis.text.y = element_text(size =12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(b.plot)
    dev.off()
}