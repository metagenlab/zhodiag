set.seed(123)
library(ggplot2)
library(tidyverse)
library(scales)

# ---------------------------
## SCRIPT ARGUMENTS ##
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
# multiqc
multiqc_data <- args[1]
# map
host_removal  <- tolower(args[2]) == "true"
host <- args[3]
# krakenuniq
run_krakenuniq <- tolower(args[4]) == "true"
krakenuniq_report <- if (run_krakenuniq) args[5] else NULL
krakenuniq_kingdom <- if (run_krakenuniq) args[6] else NULL
krakenuniq_removal <- if (run_krakenuniq) args[7] else NULL
# kraken2
run_kraken2 <- tolower(args[8]) == "true"
kraken2_report <- if (run_kraken2) args[9] else NULL
kraken2_kingdom <- if (run_kraken2) args[10] else NULL
kraken2_removal <- if (run_kraken2) args[11] else NULL
# classification by mapping
run_mapping <- tolower(args[12]) == "true"
db_name <- args[13]
bowtie2_kingdom <- if (run_mapping) args[14] else NULL
bowtie2_removal <- if (run_mapping) args[15] else NULL

# ---------------------------
## COLOURS ##
# ---------------------------
stat.colors <- c(
    "InitialReads" = "#ffffff",
    "fastp_PassedFilter" = "#f9f9f9",
    "bowtie2_HostReads" = "#969696",
    "bowtie2_nonHostReads" = "#08306b",
    "krakenuniq_ClassifiedReads" = "#ce1256",
    "krakenuniq_UnclassifiedReads" = "#b5a69a",
    "krakenuniq_HumanReads" = "#969696",
    "krakenuniq_nonHumanReads" = "#2171b5",
    "krakenuniq_Bacteria" = "#238443",
    "krakenuniq_Viruses" = "#df65b0",
    "krakenuniq_Fungi" = "#982aa8",
    "krakenuniq_removedReadsFromPlots" = "#d9d9d9",
    "krakenuniq_beforeTaxa" = "#d95f0e",
    "krakenuniq_afterTaxa" = "#ef3b2c",
    "kraken2_ClassifiedReads" = "#ce017e",
    "kraken2_UnclassifiedReads" = "#b5a69a",
    "kraken2_HumanReads" = "#969696",
    "kraken2_nonHumanReads" = "#6baed6",
    "kraken2_Bacteria" = "#addd8e",
    "kraken2_Viruses" = "#f768a1",
    "kraken2_Fungi" = "#984ea3",
    "kraken2_removedReadsFromPlots" = "#d9d9d9",
    "kraken2_beforeTaxa" = "#f2853d",
    "kraken2_afterTaxa" = "#fb6a4a",
    "bowtie2_ClassifiedReads" = "#f53b7e", 
    "bowtie2_UnclassifiedReads" = "#b5a69a",
    "bowtie2_HumanClassifiedReads" = "#969696",
    "bowtie2_nonHumanClassifiedReads" = "#1d91c0",
    "bowtie2_Bacteria" = "#78eb9e",
    "bowtie2_Viruses" = "#f582b0",
    "bowtie2_Eukaryota" = "#7b3985",
    "bowtie2_removedReadsFromPlots" = "#d9d9d9",
    "bowtie2_beforeTaxa" = "#de6b1f",
    "bowtie2_afterTaxa" = "#ed655a"
)


# ---------------------------
## TRIM STATS ##
# ---------------------------
trim_data <- read.table(
    file.path(multiqc_data, "fastp_filtered_reads_plot.txt"),
    header = TRUE, sep = '\t'
)
stats = trim_data %>% 
        column_to_rownames("Sample") %>%
        mutate(InitialReads = rowSums(across(where(is.numeric))) / 2,
                fastp_PassedFilter = Passed.Filter / 2) %>%
        select(InitialReads, fastp_PassedFilter) %>%
        rownames_to_column("Sample")

# ---------------------------
## HOST STATS ##
# ---------------------------
if (host_removal) {
    print("Mapper is bowtie2")
    map_data <- read.table(
        file.path(multiqc_data, "bowtie2_pe_plot.txt"),
        header = TRUE, sep = '\t'
    )
    map_data <- map_data[grepl(host, map_data$Sample), ] # this is to select the host map stats from this file; because if running mapping vs full db, this is added to this file too
                  
    map_data = map_data %>%
                mutate(Sample = sub(paste0("_", host), "", Sample),
                        bowtie2_HostReads = PE.mapped.uniquely,
                        bowtie2_nonHostReads = PE.mapped.discordantly.uniquely + 
                                                PE.one.mate.mapped.uniquely + 
                                                PE.one.mate.multimapped + 
                                                PE.neither.mate.aligned
                ) %>%
                select(Sample, bowtie2_HostReads, bowtie2_nonHostReads)
    stats = left_join(stats, map_data, by = "Sample")
}


# ---------------------------
## TAXONOMY
# ---------------------------

# ---------------------------
### krakenuniq
if (run_krakenuniq) {
    print("Taxonomy with krakenuniq")
    ku_data <- read.table(
        krakenuniq_report,
        header = TRUE, sep = '\t', quote = "", comment.char=""
    )
    krakenuniq = ku_data %>% group_by(sample) %>%
        summarise(krakenuniq_ClassifiedReads = sum(taxReads[taxName != "unclassified"]),
            krakenuniq_UnclassifiedReads = sum(taxReads[taxName == "unclassified"]),
            krakenuniq_nonHumanReads = sum(taxReads[!taxName %in% c("unclassified", "Homo sapiens")]),
            krakenuniq_HumanReads = sum(taxReads[taxName == "Homo sapiens"])
            ) %>%
        rename(Sample = sample)
    krakenuniq$Sample <- as.character(krakenuniq$Sample)
    stats = left_join(stats, krakenuniq, by = "Sample")

    ku_kingdoms <- read.table(
        krakenuniq_kingdom,
        header = TRUE, sep = '\t'
    )
    ku_removed <- read.table(
        krakenuniq_removal,
        header = TRUE, sep = '\t'
    )
    ku_extra = left_join(ku_kingdoms, ku_removed, by = "sample") %>%
            select(sample, krakenuniq_Bacteria, krakenuniq_Viruses, 
                    krakenuniq_Fungi, krakenuniq_removedReadsFromPlots,
                    krakenuniq_beforeTaxa, krakenuniq_afterTaxa) %>%
            rename(Sample = sample)
    ku_extra$Sample <- as.character(ku_extra$Sample)
    stats = left_join(stats, ku_extra, by = "Sample")

}
# ---------------------------
### kraken2
if (run_kraken2) {
    print("Taxonomy with kraken2")
    k2_data = read.table(
        kraken2_report,
        header = TRUE, sep = '\t', quote = "", comment.char=""
    )
    kraken2 = k2_data %>% group_by(sample) %>%
        summarise(kraken2_ClassifiedReads = sum(directCounts[taxonomy != "unclassified"]),
            kraken2_UnclassifiedReads = sum(directCounts[taxonomy == "unclassified"]),
            kraken2_nonHumanReads = sum(directCounts[!taxonomy %in% c("unclassified", "Homo sapiens")]),
            kraken2_HumanReads = sum(directCounts[taxonomy == "Homo sapiens"])
            ) %>%
        rename(Sample = sample)
    kraken2$Sample <- as.character(kraken2$Sample)
    stats = left_join(stats, kraken2, by = "Sample")

    k2_kingdoms = read.table(
        kraken2_kingdom,
        header = TRUE, sep = '\t'
    )
    k2_removed = read.table(
        kraken2_removal,
        header = TRUE, sep = '\t'
    )
    k2_extra = left_join(k2_kingdoms, k2_removed, by = 'sample') %>%
            select(sample, kraken2_Bacteria, kraken2_Viruses, 
            kraken2_Fungi, kraken2_removedReadsFromPlots,
            kraken2_beforeTaxa, kraken2_afterTaxa) %>%
            rename(Sample = sample)
    k2_extra$Sample <- as.character(k2_extra$Sample)
    stats = left_join(stats, k2_extra, by = 'Sample')

}

# ---------------------------
# Mapping bowtie2
if (run_mapping) {
    print("Taxonomy by mapping with bowtie2")
    map_data <- read.table(
        file.path(multiqc_data, "bowtie2_pe_plot.txt"),
        header = TRUE, sep = '\t'
    )
    map_data <- map_data[grepl(db_name, map_data$Sample), ] # this is to select the NON-host map stats from this file; because if running mapping vs full db, this is added to this file too
    print(db_name)
    print(map_data)
    map_data = map_data %>%
                mutate(Sample = sub(paste0("_", db_name), "", Sample),
                        bowtie2_ClassifiedReads = PE.mapped.uniquely + 
                                                    PE.mapped.discordantly.uniquely + 
                                                    PE.one.mate.mapped.uniquely + 
                                                    PE.one.mate.multimapped,
                        bowtie2_UnclassifiedReads = PE.neither.mate.aligned
                ) %>%
                select(Sample, bowtie2_ClassifiedReads, bowtie2_UnclassifiedReads)
    stats = left_join(stats, map_data, by = "Sample")
    # non-human
    db_name_nh = "filtered_noHuman"
    map_data2 <- read.table(
        file.path(multiqc_data, "samtools-flagstat-table.txt"),
        header = TRUE, sep = '\t'
    )
    map_data2 <- map_data2[grepl(db_name_nh, map_data2$Sample), ]
    map_data2 = map_data2 %>%
                mutate(Sample = sub(paste0("_", db_name_nh), "", Sample),
                        bowtie2_nonHumanClassifiedReads = Total.Reads * 1000000) %>%
                select(Sample, bowtie2_nonHumanClassifiedReads)
    stats = left_join(stats, map_data2, by = "Sample") %>%
        mutate(bowtie2_HumanClassifiedReads = bowtie2_ClassifiedReads - bowtie2_nonHumanClassifiedReads)

    b2_kingdoms <- read.table(
        bowtie2_kingdom,
        header = TRUE, sep = '\t'
    )
    b2_removed <- read.table(
        bowtie2_removal,
        header = TRUE, sep = '\t'
    )
    b2_extra <- left_join(b2_kingdoms, b2_removed, by = 'sample') %>%
            select(sample, bowtie2_Bacteria, bowtie2_Viruses, 
                    bowtie2_Eukaryota, bowtie2_removedReadsFromPlots,
                    bowtie2_beforeTaxa, bowtie2_afterTaxa) %>%
            rename(Sample = sample)
    b2_extra$Sample <- as.character(b2_extra$Sample)
    stats = left_join(stats, b2_extra, by = "Sample")
}

# ---------------------------
## SAVE REPORT TABLE
# ---------------------------
write.table(stats, 
            'stats_report.tsv', 
            quote = FALSE, sep = '\t', 
            col.names = TRUE, row.names = FALSE)


# ---------------------------
## PLOT STATS
# ---------------------------
stat.plot <- stats %>% pivot_longer(cols = -Sample)
stat.plot$name <- factor(stat.plot$name, levels = colnames(stats))
p = ggplot(stat.plot, aes(x = name, y = log10(value), fill = name, label = value)) +
  geom_bar(stat = 'identity', colour = "#d9d9d9", size = 0.005) +
  geom_text(angle = 90, size = 2, hjust = -0.2, colour = '#000000') + 
  facet_wrap(facets = "Sample", ncol = 4) +
  scale_fill_manual(values = stat.colors) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = '', y = 'Number of reads (log10)', fill = '') +
  theme(axis.text.x = element_blank(),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))

# n_samples = length(colnames(stats))
# plot_width = n_samples * 4
pdf("stats_report.pdf", width = 14)
print(p)
dev.off()
