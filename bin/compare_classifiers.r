set.seed(123)
library(ggplot2)
library(tidyverse)
library(scales)
library(stringr)
sessionInfo()

# ---------------------------
## SCRIPT ARGUMENTS ##
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

run_krakenuniq <- tolower(args[1]) == "true"
krakenuniq_df <- args[2]

run_kraken2 <- tolower(args[3]) == "true"
kraken2_df <- args[4]

run_bowtie2 <- tolower(args[5]) == "true"
bowtie2_df <- args[6]

# ---------------------------
## READ DATA ##
# ---------------------------
full_df <- data.frame()

if (run_krakenuniq) {
    ku <- read.table( krakenuniq_df,
        header = TRUE, sep = '\t'
    )
    ku.l <- ku %>% 
        rename(taxid = taxID,
                species = taxName) %>%
        pivot_longer(
        cols = -c(taxid, species),
        names_to = "Sample",
        values_to = "Reads"
    ) %>% 
    mutate(classifier = 'KrakenUniq')

    full_df <- rbind(full_df, ku.l) 
}

if (run_kraken2) {
    k2 <- read.table( kraken2_df,
        header = TRUE, sep = '\t'
    )
    k2.l <- k2 %>% 
        rename(species = taxonomy) %>%
        pivot_longer(
        cols = -c(taxid, species),
        names_to = "Sample",
        values_to = "Reads"
    ) %>% 
    mutate(classifier = 'Kraken2')

    full_df <- rbind(full_df, k2.l) 
}

if (run_bowtie2) {
    b2 <- read.table( bowtie2_df,
        header = TRUE, sep = '\t'
    )
    b2.l <- b2 %>% pivot_longer(
        cols = -c(taxid, species),
        names_to = "Sample",
        values_to = "Reads"
    ) %>% 
    mutate(classifier = 'Bowtie2')

    full_df <- rbind(full_df, b2.l) 
}

write.table(full_df, "classifiers_combined.tsv",
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)

# ---------------------------
## PER SAMPLE ##
# ---------------------------
# plot size
height_per_species <- 0.2  # inches per species
base_height <- 4           # minimal height in inches
width_per_sample <- 0.5    # inches per sample
base_width <- 8            # minimal width in inches

for (s in unique(full_df$Sample)) {

    # Filter for current sample and species with at least one read
    df.s <- full_df %>% 
        filter(Sample == s) %>%
        group_by(taxid, species) %>%
        filter(sum(Reads) > 0) %>%
        ungroup()

    # Compute total reads per species
    totals <- df.s %>%
        group_by(taxid, species) %>%
        summarise(totReads_species = sum(Reads), .groups = "drop")

    # Join totals to df.s
    df.s <- df.s %>%
        left_join(totals, by = c("taxid", "species"))

    # Pivot to wide table
    df.w <- df.s %>%
        pivot_wider(
            id_cols = c(taxid, species, totReads_species),
            names_from = classifier,
            values_from = Reads,
            values_fill = 0
        ) %>%
        arrange(desc(totReads_species)) %>%
        select(-totReads_species)

    # Write wide table
    write.table(df.w,
                paste0(s, "_classifier_comparison.tsv"),
                col.names = TRUE,
                row.names = FALSE,
                sep = '\t',
                quote = FALSE)

    # Prepare plot size
    n_species <- nrow(df.w)
    n_samples <- length(unique(df.s$classifier))
    plot_height <- max(base_height, n_species * height_per_species)
    plot_width  <- max(base_width,  n_samples * width_per_sample)

    # Heatmap plot (species sorted by total reads)
    p <- ggplot(df.s,
                aes(x = classifier,
                    y = reorder(species, totReads_species),
                    fill = Reads,
                    label = Reads)) +
        geom_tile() +
        geom_text(colour = "white") +
        labs(x = '', y = '') +
        ggtitle(s) +
        theme_classic()

    pdf(paste0(s, "_classifier_comparison.pdf"),
        height = plot_height, width = plot_width)
    print(p)
    dev.off()
}