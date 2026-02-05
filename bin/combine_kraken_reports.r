set.seed(123)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
meta <- args[1]
outfile_prefix <- args[2]
file_list <- args[3:length(args)]

metadata <- fread(meta)
metadata = metadata %>% select(sample, group)
metadata$sample = as.character(metadata$sample)

df = data.frame()

for(f in file_list){
    s = sub("_db.*", "", f)
    a = fread(f, header = FALSE, sep = '\t', strip.white = TRUE, fill = TRUE)
    b = a %>% mutate(sample = s)
    df = bind_rows(df, b)
}
df$sample = as.character(df$sample)
dfg = left_join(df, metadata, by = "sample")
colnames(dfg) = c("perc_reads", "totalCounts", "directCounts", 
                    "minimizers", "distinctMinimizers", "rank",
                    "taxid", "taxonomy", "sample", "group")
dfg = dfg %>% mutate(taxonomy = gsub("#", "", taxonomy))

write.table(dfg, paste0(outfile_prefix),
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)