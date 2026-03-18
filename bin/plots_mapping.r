#!/usr/bin/env Rscript
set.seed(123)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
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
min_reads <- as.numeric(args[3])
contaminants <- args[4]
file_list <- args[5:length(args)]


# ---------------------------
#  METADATA
# ---------------------------
metadata <- fread(meta)
metadata = metadata %>% select(sample, group)
metadata$sample = as.character(metadata$sample)

# ---------------------------
# CONCATENATE SAMPLE TABLES; ADD METADATA
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

write.table(dfg, paste0('bowtie2_combined_table_by_', level, '.tsv'),
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)

# ---------------------------
# PLOTS FROM TAXID TABLE
# ---------------------------
cov_colors <- c(
  "0-300"      = "#f7fbff",
  "300-500"    = "#deebf7",
  "500-1000"   = "#c6dbef",
  "1000-2000"  = "#9ecae1",
  "2000-5000"  = "#6baed6",
  "5000-10000" = "#3182bd",
  "10000-100000" = "#08519c",
  "100000+"    = "#0c2c84"
)
height_per_species <- 0.2  # inches per species
base_height <- 10           # minimal height in inches
width_per_sample <- 1    # inches per sample
base_width <- 10            # minimal width in inches

if(level == 'taxid') {
  print(paste0("Taxa with less than ", min_reads, " reads were removed from plots"))

  if (!is.null(contaminants)) {
    contams <- read.table(contaminants, header = FALSE, sep = '\t')
    colnames(contams) <- c('contaminant', 'code')
    print("The following taxa were removed: ")
    print(contams)
  }

  # WIDE TABLE OF READS PER SAMPLE AND SPECIES
  w.dfg <- dfg %>%
    select(taxid, mappedReads, species, sample, group) %>%
    unite(sample_group, sample, group, sep = "_") %>%
    pivot_wider(
      names_from = sample_group,
      values_from = mappedReads,
      values_fill = 0
    )
  write.table(w.dfg, paste0("bowtie2_read_table_at_species_level.tsv"),
    row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


  # SUMMARY BY KINGDOM
  expected_domains <- c(
    "bowtie2_Viruses",
    "bowtie2_Bacteria",
    "bowtie2_Eukaryota"
  )
  sum.domains <- dfg %>% 
    mutate(domain = case_when(
      domain == "Viruses" ~ "bowtie2_Viruses",
      domain == "Bacteria" ~ "bowtie2_Bacteria",
      domain == "Eukaryota" ~ "bowtie2_Eukaryota",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(domain)) %>%
    mutate(domain = factor(domain, levels = expected_domains)) %>%
    group_by(sample, domain) %>%
    summarise(reads = sum(mappedReads), .groups = "drop") %>%
    complete(sample, domain = expected_domains, fill = list(reads = 0)) %>%
    pivot_wider(names_from = domain, values_from = reads, values_fill = 0)
  write.table(sum.domains, "bowtie2_summary_kingdoms.tsv",
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')

  # FILTERS:
  #        * have at least "min_reads" reads
  #        * exclude taxa specified in params.contaminants file
  #        * only detected in control group
  bef.reads <- dfg %>%
    group_by(sample) %>%
    summarise(beforeReads = sum(mappedReads),
              bowtie2_beforeTaxa = n())
  
  dfg.filt <- dfg %>%
    select(taxid, nBases_covered, mappedReads, mean_cov, species, sample, group) %>%
    filter(mappedReads > min_reads)
  
  if (!is.null(contaminants)) {
    sp.filter = contams %>% filter(code == 'species') %>% pull(contaminant)
    gn.filter = contams %>% filter(code == 'genus') %>% pull(contaminant)
    tax.filter = contams %>% filter(code == 'taxid') %>% pull(contaminant)

    # save removed contaminants
    dfg.removed <- dfg.filt %>% 
        filter(species %in% sp.filter |
              Reduce('|', lapply(gn.filter, function(g) startsWith(species, g))) |
              taxid %in% tax.filter)
    write.table(dfg.removed, "bowtie2_removedContaminants.tsv",
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    print(dfg.removed)

    # remove the contaminants
    dfg.filt <- dfg.filt %>% 
        filter(!species %in% sp.filter,
                !Reduce('|', lapply(gn.filter, function(g) startsWith(species, g))),
                !taxid %in% tax.filter
                )
      }

  a <- dfg.filt %>%
    group_by(species) %>%
    mutate(
      control_reads = sum(mappedReads[group == "control"]),
      other_reads = sum(mappedReads[group != "control"])
    ) %>% 
    ungroup() %>%
    filter(!(control_reads > 0 & other_reads == 0)) %>%
    select(-control_reads, -other_reads)
  
  aft.reads <- a %>%
  group_by(sample) %>%
  summarise(afterReads = sum(mappedReads),
            bowtie2_afterTaxa = n())
  
  rem.reads <- left_join(bef.reads, aft.reads, by = "sample") %>%
    mutate(bowtie2_removedReadsFromPlots = beforeReads - afterReads)
  write.table(rem.reads, "bowtie2_removedReadsFromPlots.tsv",
    row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  
  #  WIDE TABLE OF READS PER SAMPLE AND SPECIES POST FILTERING
  watot <- a %>%
    select(taxid, mappedReads, species, sample, group) %>%
    unite(sample_group, sample, group, sep = "_") %>%
    pivot_wider(
      names_from = sample_group,
      values_from = mappedReads, 
      values_fill = 0
    )
  write.table(watot, "bowtie2_read_table_at_species_level_postCleaning.tsv",
        row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  

  # REORDER BY TOTAL READS PER SPECIES ACROSS ALL SAMPLES
  totals <- a %>%
    group_by(species) %>%
    summarise(totReads_species = sum(mappedReads))
  a <- a %>% left_join(totals, by = 'species')
  a$species <- factor(a$species, levels = totals$species[order(totals$totReads_species)])


  # add coverage ranges for plotting  colours
  a <- a %>%
    mutate(coverage = cut(
      nBases_covered,
      breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, 100000, Inf),
      labels = c("0-300", "300-500", "500-1000", "1000-2000",
                  "2000-5000", "5000-10000", "10000-100000", "100000+"),
      right = FALSE
    ))
  # plot size
  n_species <- length(unique(a$species))
  n_samples <- length(unique(a$sample))
  plot_height <- max(base_height, n_species * height_per_species)
  plot_width  <- max(base_width,  n_samples * width_per_sample)

  # HETMAP SPECIES: BASES COVERED and MAPPED READS
  pdf('heatmap_reads_fillByCoverage.pdf', height = plot_height, width = plot_width)
  p = ggplot(a, 
    aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
    geom_tile() +
    geom_text(colour='black') +
    scale_fill_manual(values = cov_colors) +
    labs(x = '', y = '') +
    facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()

  # HEATMAPS BY GROUP
  for(gr in setdiff(unique(a$group), unique(a$group[grepl("^control", a$group)]))){
    print(gr)
    control_groups <- unique(a$group[grepl("^control", a$group)])
    dtpm <- a %>% filter(group %in% c(gr, control_groups))

    # filter out taxa only in control vs group
    dtp <- dtpm %>%
    group_by(species) %>%
    mutate(
      control_reads = sum(mappedReads[group == "control"]),
      other_reads   = sum(mappedReads[group != "control"])
    ) %>%
    ungroup() %>%
    filter(!(control_reads > 0 & other_reads == 0)) %>%
    select(-control_reads, -other_reads)

    # reorder by total reads per species across all samples
    totals.gr <- dtp %>%
      group_by(species) %>%
      summarise(totReads_species = sum(mappedReads))
    dtp <- dtp %>% left_join(totals.gr, by = "species")
    dtp$species <- factor(dtp$species, levels = totals.gr$species[order(totals.gr$totReads_species)])

    # plot size
    n_species <- length(unique(dtp$species))
    n_samples <- length(unique(dtp$sample))
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    pdf(paste0(gr, "_group_heatmap_reads_fillByCoverage.pdf"), height = plot_height, width = plot_width)
    p = ggplot(dtp, 
      aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
      geom_tile() +
      geom_text(colour='black') +
      scale_fill_manual(values = cov_colors) +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()

  # HEATMAPS BY SAMPLE
  for(gr in setdiff(unique(a$sample), unique(a$sample[grepl("^control", a$group)]))){
    print(gr)
    control_groups <- unique(a$sample[grepl("^control", a$group)])
    dtpm <- a %>% filter(sample %in% c(gr, control_groups))
    
    # filter out taxa only in control vs sample
    dtp <- dtpm %>%
      group_by(species) %>%
      mutate(
        control_reads = sum(mappedReads[group == "control"]),
        other_reads   = sum(mappedReads[group != "control"])
      ) %>%
      ungroup() %>%
      filter(!(control_reads > 0 & other_reads == 0)) %>%
      select(-control_reads, -other_reads)
    
    # reorder by total reads per species across all samples
    totals.gr <- dtp %>%
      group_by(species) %>%
      summarise(totReads_species = sum(mappedReads))
    dtp <- dtp %>% left_join(totals.gr, by = "species")
    dtp$species <- factor(dtp$species, levels = totals.gr$species[order(totals.gr$totReads_species)])
    
    # plot size
    n_species <- length(unique(dtp$species))
    n_samples <- length(unique(dtp$sample))
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)
    
    pdf(paste0(gr, "_sample_heatmap_reads_fillByCoverage.pdf"), height = plot_height, width = plot_width)
    p = ggplot(dtp, 
               aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
      geom_tile() +
      geom_text(colour='black') +
      scale_fill_manual(values = cov_colors) +
      labs(x = '', y = '') +
      facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)
    dev.off()
  }
}

  # READS VS COVERED POSITIONS FOR EACH SAMPLE
  for(s in unique(a$sample)){
    print(s)
    pdf(paste0(s, "_reads_vs_basesCovered.pdf"), 
    width = 10, height = 10)
    p = ggplot(a %>% filter(sample == s), 
          aes(x = log2(mappedReads), y = log2(nBases_covered), label = species, size = mean_cov)) +
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
}

  # # VIRUS
  #   virus <- dfg %>% filter(domain == "Viruses") %>%
  #     group_by(sample, species, group) %>%
  #     summarise(
  #       mappedReads = sum(mappedReads, na.rm = TRUE),
  #       nBases_covered = sum(nBases_covered, na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     filter(!is.na(species)) %>%
  #     # filter(mappedReads > 1, nBases_covered > 151) %>%
  #     group_by(species) %>%
  #     mutate(total_reads = sum(mappedReads, na.rm = TRUE)) %>%
  #     ungroup() %>%
  #     mutate(species = reorder(species, total_reads)) %>%
  #     mutate(coverage = cut(
  #       nBases_covered,
  #       breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
  #       labels = c("0-300", "300-500", "500-1000", "1000-2000",
  #                 "2000-5000", "5000-10000", "10000+"),
  #       right = FALSE
  #     ))
      
  #   # plot size
  #   n_species <- length(unique(virus$species))
  #   n_samples <- length(unique(virus$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)


  #   pdf(paste0('heatmap_VIRUSES_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
  #   v.plot = ggplot(virus , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour = 'black') +
  #   facet_grid(~group, scales = 'free_x', space = 'free') +
  #   labs(x = '', y = '') +
  #   theme_classic() +
  #   scale_fill_manual(values = cov_colors) +
  #   theme(axis.text.y = element_text(size =12),
  #         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(v.plot)
  #   dev.off()

  #   # filter taxids detected only in control
  #   df.filt.v <- virus  %>%
  #     group_by(species) %>% 
  #     filter(
  #       any(mappedReads > 0 & group != "control")
  #     ) %>%
  #     ungroup()

  #   # plot size
  #   n_species <- length(unique(df.filt.v$species))
  #   n_samples <- length(unique(df.filt.v$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)

  #   # HETMAP SPECIES: BASES COVERED and MAPPED READS
  #   pdf(paste0('heatmap_VIRUSES_fillCoverage_labelMappedReads_filtered.pdf'), height = plot_height, width = plot_width)
  #   p = ggplot(df.filt.v, 
  #   aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour='black') +
  #   scale_fill_manual(values = cov_colors) +
  #   labs(x = '', y = '') +
  #   facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(p)
  #   dev.off()


  # # EUKARYOTA
  #   euka <-  dfg %>% filter(domain == "Eukaryota") %>%
  #     group_by(sample, species, group) %>%
  #     summarise(
  #       mappedReads = sum(mappedReads, na.rm = TRUE),
  #       nBases_covered = sum(nBases_covered, na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     filter(!is.na(species)) %>%
  #     # filter(mappedReads > 2, nBases_covered > 151) %>%
  #     group_by(species) %>%
  #     mutate(total_reads = sum(mappedReads, na.rm = TRUE)) %>%
  #     ungroup() %>%
  #     mutate(species = reorder(species, total_reads)) %>%
  #     mutate(coverage = cut(
  #       nBases_covered,
  #       breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
  #       labels = c("0-300", "300-500", "500-1000", "1000-2000",
  #                 "2000-5000", "5000-10000", "10000+"),
  #       right = FALSE
  #     ))
    
  #   # plot size
  #   n_species <- length(unique(euka$species))
  #   n_samples <- length(unique(euka$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)


  #   pdf(paste0('heatmap_EUKARYOTA_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
  #   e.plot = ggplot(euka , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour = 'black') +
  #   facet_grid(~group, scales = 'free_x', space = 'free') +
  #   labs(x = '', y = '') +
  #   scale_fill_manual(values = cov_colors) +
  #   theme_classic() +
  #   theme(axis.text.y = element_text(size =12),
  #         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(e.plot)
  #   dev.off()

  #   # filter taxids detected only in control
  #   df.filt.e <- euka  %>%
  #     group_by(species) %>% 
  #     filter(
  #       any(mappedReads > 0 & group != "control")
  #     ) %>%
  #     ungroup()

  #   # plot size
  #   n_species <- length(unique(df.filt.e$species))
  #   n_samples <- length(unique(df.filt.e$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)

  #   # HETMAP SPECIES: BASES COVERED and MAPPED READS
  #   pdf(paste0('heatmap_EUKARYOTA_fillCoverage_labelMappedReads_filtered.pdf'), height = plot_height, width = plot_width)
  #   p = ggplot(df.filt.e, 
  #   aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour='black') +
  #   scale_fill_manual(values = cov_colors) +
  #   labs(x = '', y = '') +
  #   facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(p)
  #   dev.off()



  # # BACTERIA
  #   bact <-   dfg %>% filter(domain == "Bacteria") %>%
  #     group_by(sample, species, group) %>%
  #     summarise(
  #       mappedReads = sum(mappedReads, na.rm = TRUE),
  #       nBases_covered = sum(nBases_covered, na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     filter(!is.na(species)) %>%
  #     # filter(mappedReads > 2, nBases_covered > 151) %>%
  #     group_by(species) %>%
  #     mutate(total_reads = sum(mappedReads, na.rm = TRUE)) %>%
  #     ungroup() %>%
  #     mutate(species = reorder(species, total_reads)) %>%
  #     mutate(coverage = cut(
  #       nBases_covered,
  #       breaks = c(0, 300, 500, 1000, 2000, 5000, 10000, Inf),
  #       labels = c("0-300", "300-500", "500-1000", "1000-2000",
  #                 "2000-5000", "5000-10000", "10000+"),
  #       right = FALSE
  #     ))

  #   # plot size
  #   n_species <- length(unique(bact$species))
  #   n_samples <- length(unique(bact$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)

  #   pdf(paste0('heatmap_BACTERIA_fillCoverage_labelMappedReads.pdf'), height = plot_height, width = plot_width)
  #   b.plot = ggplot(bact , aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour = 'black') +
  #   facet_grid(~group, scales = 'free_x', space = 'free') +
  #   labs(x = '', y = '') +
  #   scale_fill_manual(values = cov_colors) +
  #   theme_classic() +
  #   theme(axis.text.y = element_text(size =12),
  #         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(b.plot)
  #   dev.off()

  #   # filter taxids detected only in control
  #   df.filt.b <- bact  %>%
  #     group_by(species) %>% 
  #     filter(
  #       any(mappedReads > 0 & group != "control")
  #     ) %>%
  #     ungroup()

  #   # plot size
  #   n_species <- length(unique(df.filt.b$species))
  #   n_samples <- length(unique(df.filt.b$sample))
  #   plot_height <- max(base_height, n_species * height_per_species)
  #   plot_width  <- max(base_width,  n_samples * width_per_sample)

  #   # HETMAP SPECIES: BASES COVERED and MAPPED READS
  #   pdf(paste0('heatmap_BACTERIA_fillCoverage_labelMappedReads_filtered.pdf'), height = plot_height, width = plot_width)
  #   p = ggplot(df.filt.b, 
  #   aes(x = sample, y = species, fill = coverage, label = mappedReads)) +
  #   geom_tile() +
  #   geom_text(colour='black') +
  #   scale_fill_manual(values = cov_colors) +
  #   labs(x = '', y = '') +
  #   facet_grid(.~factor(group), scales = 'free_x', space = 'free') +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #   print(p)
  #   dev.off()


