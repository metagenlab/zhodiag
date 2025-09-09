#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)

# ---------------------------
# Parse command-line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
map_summary <- args[2]
output_prefix <- args[3]

# ---------------------------
# Load input
# ---------------------------
print('Read depth file')
df1 <- data.table::fread(filename, header = FALSE, sep = '\t')
colnames(df1) <- c("id", "pos", "depth")

# Split id into accession and taxid
df1 <- df1 %>%
  separate(col = id, into = c("accession", "taxid"),
           sep = "\\|", remove = TRUE)

# ---------------------------
# Load mapping summary
# ---------------------------
print("Read map summary file")
df2 <- data.table::fread(map_summary, header = FALSE, sep = '\t')
colnames(df2) <- c("nReads", "id")

# Split id into accession and taxid
df2 <- df2 %>%
  separate(col = id, into = c("accession", "taxid"),
           sep = "\\|", remove = TRUE)

# ---------------------------
# Summaries by accession
# ---------------------------
print("Table by accession")
summary_by_accession <- df1 %>%
  group_by(accession) %>%
  summarise(
    taxid = first(taxid),
    length = max(pos),
    fraction_covered = sum(depth > 0) / max(pos),
    mean_depth = if (any(depth > 0)) mean(depth[depth > 0]) else 0,
    median_depth = if (any(depth > 0)) median(depth[depth > 0]) else 0,
    .groups = "drop"
  ) %>%
  arrange(desc(fraction_covered))

# add map summary
summary_by_accession <- left_join(summary_by_accession, df2, by = c('accession', 'taxid'))

write.table(summary_by_accession,
            paste0(output_prefix, "summary_statistics_by_accession.tsv"),
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)


# ---------------------------
# Summaries by taxid
# ---------------------------
print("Table by taxID")
summary_by_taxid <- summary_by_accession %>%
  group_by(taxid) %>%
  summarise(
    nContigs = n(),
    length = sum(length),
    mean_cov = mean(fraction_covered),
    mean_depth = mean(mean_depth),
    median_depth = median(median_depth),
    totalReads = sum(nReads),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cov))

write.table(summary_by_taxid,
            paste0(output_prefix, "summary_statistics_by_taxid.tsv"),
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)
