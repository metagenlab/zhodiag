set.seed(123)
library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
# filename <- args[1]
multiqc_data <- args[1]
trim_tool <- args[2]
mapper <- args[3]
run_krakenuniq <- tolower(args[4]) == "true"
krakenuniq_report <- if (run_krakenuniq) args[5] else NULL
run_kraken2 <- tolower(args[6]) == "true"
kraken2_report <- if (run_kraken2) args[7] else NULL

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

## HOST STATS ##
    print("Mapper is bowtie2")
    if (mapper == "bowtie2") {
        map_data <- read.table(
            file.path(multiqc_data, "bowtie2_pe_plot.txt"),
            header = TRUE, sep = '\t'
        )
    map_data = map_data %>%
                mutate(Sample = sub("_[^_]*$", "", Sample),
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
## SAVE REPORT TABLE
    write.table(stats, 
                'stats_report.tsv', 
                quote = FALSE, sep = '\t', 
                col.names = TRUE, row.names = FALSE)
} 


## PLOT STATS
stat.plot <- stats %>% pivot_longer(cols = -Sample)
stat.plot$name <- factor(stat.plot$name, levels = colnames(stats))
p = ggplot(stat.plot, aes(x = name, y = value, fill = name)) +
  geom_bar(stat = 'identity') +
  facet_wrap(facets = "Sample") +
  theme_bw() +
  labs(x = '', y = 'Number of reads', fill = '') +
  theme(axis.text.x = element_blank()) +
  scale_fill_viridis_d()

# n_samples = length(colnames(stats))
# plot_width = n_samples * 4
pdf("stats_report.pdf")
print(p)
dev.off()