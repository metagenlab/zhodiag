#!/usr/bin/env Rscript
set.seed(123)
suppressPackageStartupMessages({
  library(tidyverse)
  library(taxonomizr)
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
filename <- args[1]
map_summary <- args[2]
output_prefix <- args[3]

# ---------------------------
# Path to taxonomy DB
# ---------------------------
dbDir  <- file.path(getwd(), "data")
dbFile <- file.path(dbDir, "accessionTaxa.sql")

dbFile <- Sys.getenv("TAXONOMY_DB", dbFile)

# Ensure data directory exists
if (!dir.exists(dbDir)) {
  dir.create(dbDir, recursive = TRUE)
}

# ---------------------------
# Build taxonomy DB if needed
# ---------------------------
if (!file.exists(dbFile)) {
  message("Taxonomy database not found at: ", dbFile)
  message("Building it now... (this may take several minutes)")
  prepareDatabase(dbFile)
} else {
  message("Using existing taxonomy database: ", dbFile)
}

# ---------------------------
# Load input
# ---------------------------
print('Read depth file')
# df1 <- data.table::fread(filename, header = FALSE, sep = '\t')
# colnames(df1) <- c("id", "pos", "depth")

# # Split id into accession and taxid
# df1 <- df1 %>%
#   separate(col = id, into = c("accession", "taxid"),
#            sep = "\\|", remove = TRUE)
summary_by_accession <- data.table::fread(filename, header = TRUE, sep = '\t')
summary_by_accession$taxid <- as.integer(summary_by_accession$taxid)
# ---------------------------
# Load mapping summary
# ---------------------------
print("Read map summary file")
df2 <- data.table::fread(map_summary, header = FALSE, sep = ' ')
colnames(df2) <- c("mappedReads", "id")

# Split id into accession and taxid
df2 <- df2 %>%
  separate(col = id, into = c("accession", "taxid"),
           sep = "\\|", remove = TRUE)
df2$taxid = as.integer(df2$taxid)
# ---------------------------
# Summaries by accession
# ---------------------------
print("Table by accession")
# summary_by_accession <- df1 %>%
#   group_by(accession) %>%
#   summarise(
#     taxid = first(taxid),
#     length = max(pos),
#     nBases_covered = sum(depth > 0),
#     fraction_covered = sum(depth > 0) / max(pos),
#     mean_depth = if (any(depth > 0)) mean(depth[depth > 0]) else 0,
#     median_depth = if (any(depth > 0)) median(depth[depth > 0]) else 0,
#     .groups = "drop"
#   ) %>%
#   arrange(desc(nBases_covered))

# add map summary
summary_by_accession <- left_join(summary_by_accession, df2, by = c('accession', 'taxid'))

taxonomy_acc <- getTaxonomy(as.integer(summary_by_accession$taxid), dbFile)
summary_by_accession <- cbind(summary_by_accession, taxonomy_acc)

write.table(summary_by_accession,
            paste0(output_prefix, "_summary_statistics_by_accession.tsv"),
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)


# # PLOT NREADS VS COVERAGE
# pdf(paste0(output_prefix, "_nMappedReads_vs_coverageFraction_byAccession.pdf"), height = 10, width = 10)
# p = ggplot(summary_by_accession, 
#             aes(x = log2(mappedReads), y = fraction_covered, colour = mean_depth, label = species)) +
#     geom_point(size = 2.5) +
#     geom_text_repel() +
#     scale_y_continuous(limits = c(0, 1)) +
#     scale_x_continuous(limits = c(0, NA)) +
#     theme_classic() +
#     labs(x = "log2(Mapped Reads)", y = "Fraction covered", colour = "mean depth") +
#     theme(axis.title = element_text(size = 14),
#           axis.text = element_text(size = 14)) +
#     ggtitle(output_prefix)    
# print(p)
# dev.off()

# PLOT NREADS VS COVERED POSITIONS
pdf(paste0(output_prefix, "_nMappedReads_vs_nBasesCovered_byAccession.pdf"), height = 10, width = 10)
p = ggplot(summary_by_accession, 
            aes(x = log2(mappedReads), y = nBases_covered, colour = fraction_covered, label = species)) +
    geom_point(size = 2.5) +
    geom_text_repel() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_classic() +
    labs(x = "log2(Mapped Reads)", y = "Bases covered", colour = "Fraction covered") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ggtitle(output_prefix) 
print(p)
dev.off()

# ---------------------------
# ---------------------------
# ---------------------------
# Summaries by taxid
# ---------------------------
# ---------------------------
# ---------------------------
summary_by_taxid <- summary_by_accession %>%
  group_by(taxid) %>%
  summarise(
    nContigs = n(),
    length = sum(length),
    nBases_covered = sum(nBases_covered),
    mean_cov = mean(fraction_covered),
    mean_depth = mean(mean_depth),
    # median_depth = median(median_depth),
    mappedReads = sum(mappedReads),
    .groups = "drop"
  ) %>%
  arrange(desc(nBases_covered))

taxonomy_tax <- getTaxonomy(as.integer(summary_by_taxid$taxid), dbFile)
summary_by_taxid <- cbind(summary_by_taxid, taxonomy_tax)

write.table(summary_by_taxid,
            paste0(output_prefix, "_summary_statistics_by_taxid.tsv"),
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)

# PLOT NREADS VS COVERED POSITIONS
pdf(paste0(output_prefix, "_nMappedReads_vs_nBasesCovered_byTaxid.pdf"), height = 10, width = 10)
p = ggplot(summary_by_taxid, 
            aes(x = log2(mappedReads), y = nBases_covered, colour = mean_cov, label = species)) +
    geom_point(size = 2.5) +
    geom_text_repel() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_classic() +
    labs(x = "log2(Mapped Reads)", y = "Bases covered", colour = "Mean Fraction\ncovered") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ggtitle(output_prefix) 
print(p)
dev.off()

# GROUP BY GENUS
summary_by_taxid_by_genus = summary_by_taxid %>% group_by(genus) %>% 
  summarise(mappedReads = sum(mappedReads),
            nBases_covered = sum(nBases_covered),
            genus = first(genus))

pdf(paste0(output_prefix, "_nMappedReads_vs_nBasesCovered_byTaxid_genusLevel.pdf"), height = 10, width = 10)
p = ggplot(summary_by_taxid_by_genus, 
            aes(x = log2(mappedReads), y = nBases_covered, label = genus)) +
    geom_point(size = 2.5) +
    geom_text_repel() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_classic() +
    labs(x = "log2(Mapped Reads)", y = "Bases covered") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ggtitle(output_prefix) 
print(p)
dev.off()
