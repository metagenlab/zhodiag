set.seed(123)
library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
# filename <- args[1]
multiqc_data <- args[1]
trim_tool <- args[2]
mapper <- args[3]
host <- args[4]
run_krakenuniq <- tolower(args[5]) == "true"
krakenuniq_report <- if (run_krakenuniq) args[6] else NULL
run_kraken2 <- tolower(args[7]) == "true"
kraken2_report <- if (run_kraken2) args[8] else NULL
run_mapClassified <- tolower(args[9]) == "true"

if (trim_tool == "fastp") {
    print("Trim tool is fastp")
## TRIM STATS ##
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
}

## HOST STATS ##
if (mapper == "bowtie2") {
    print("Mapper is bowtie2")
    map_data <- read.table(
        file.path(multiqc_data, "bowtie2_pe_plot.txt"),
        header = TRUE, sep = '\t'
    )
    map_data <- map_data[grepl(host, map_data$Sample), ] # this is to select the host map stats from this file; because if running mapping vs full db, this is added to this file too
    map_data = map_data %>%
                mutate(Sample = stats$Sample, # helps clean the samplename column, but not the best practice,,,
                        bowtie2_HostReads = PE.mapped.uniquely,
                        bowtie2_nonHostReads = PE.mapped.discordantly.uniquely + PE.one.mate.mapped.uniquely + PE.one.mate.multimapped + PE.neither.mate.aligned
                ) %>%
                select(Sample, bowtie2_HostReads, bowtie2_nonHostReads)
    stats = left_join(stats, map_data, by = "Sample")
}

## TAXONOMY
### krakenuniq
if (run_krakenuniq) {
    print("Taxonomy with krakenuniq")
    ku_data <- read.table(
        krakenuniq_report,
        header = TRUE, sep = '\t'
    )
    krakenuniq = ku_data %>% group_by(sample) %>%
        summarise(krakenuniq_ClassifiedReads = sum(taxReads[taxName != "unclassified"]),
            krakenuniq_UnclassifiedReads = sum(taxReads[taxName == "unclassified"]),
            krakenuniq_HumanReads = sum(taxReads[taxName == "Homo sapiens"]),
            krakenuniq_nonHumanReads = sum(taxReads[!taxName %in% c("unclassified", "Homo sapiens")])
            ) %>%
        rename(Sample = sample)
    stats = left_join(stats, krakenuniq, by = "Sample")

}
### kraken2
if (run_kraken2) {
    print("Taxonomy with kraken2")
    k2_data = read.table(
        kraken2_report,
        header = TRUE, sep = '\t'
    )
    kraken2 = k2_data %>% group_by(sample) %>%
        summarise(kraken2_ClassifiedReads = sum(directCounts[taxonomy != "unclassified"]),
            kraken2_UnclassifiedReads = sum(directCounts[taxonomy == "unclassified"]),
            kraken2_HumanReads = sum(directCounts[taxonomy == "Homo sapiens"]),
            kraken2_nonHumanReads = sum(directCounts[!taxonomy %in% c("unclassified", "Homo sapiens")])
            ) %>%
        rename(Sample = sample)
    stats = left_join(stats, kraken2, by = "Sample")

}
### mapping classified reads to fulldb
if (run_mapClassified) {
    print("Mapping of classified reads to full db")
    if (mapper == "bowtie2") {
        print("Mapper is bowtie2")
        map_data <- read.table(
            file.path(multiqc_data, "samtools-flagstat-table.txt"),
            header = TRUE, sep = '\t'
        )
        map_data = map_data %>%
                    mutate(Sample = stats$Sample,
                            bowtie2_nonHumanFullDbMappedReads = Total.Reads * 1000000) %>%
                    select(Sample, bowtie2_nonHumanFullDbMappedReads)
        stats = left_join(stats, map_data, by = "Sample")

    }
}



## SAVE REPORT TABLE
write.table(stats, 
            'stats_report.tsv', 
            quote = FALSE, sep = '\t', 
            col.names = TRUE, row.names = FALSE)


## PLOT STATS
stat.plot <- stats %>% pivot_longer(cols = -Sample)
stat.plot$name <- factor(stat.plot$name, levels = colnames(stats))
p = ggplot(stat.plot, aes(x = name, y = log10(value), fill = name, label = value)) +
  geom_bar(stat = 'identity') +
  geom_text(angle = 90, size = 2, hjust = 1, colour = "#666666") + 
  facet_wrap(facets = "Sample") +
  theme_classic() +
  labs(x = '', y = 'Number of reads (log10)', fill = '') +
  theme(axis.text.x = element_blank(),
        legend.text  = element_text(size = 7)) +
  scale_fill_viridis_d()

# n_samples = length(colnames(stats))
# plot_width = n_samples * 4
pdf("stats_report.pdf")
print(p)
dev.off()