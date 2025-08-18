library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
meta <- args[1]
outfile_prefix <- args[2]
file_list <- args[3:length(args)]

metadata = read.table(meta, header = TRUE, sep =',')
metadata = metadata %>% select(sample, group)
metadata$sample = as.character(metadata$sample)

df = data.frame()

for(f in file_list){
    s = sub("_.*", "", f)
    a = read.table(f, header = FALSE, sep = '\t', strip.white = TRUE, fill = TRUE)
    a = a %>% mutate(sample = s)
    df = rbind(df, a)
}
df$sample = as.character(df$sample)
dfg = left_join(df, metadata, by = "sample")
colnames(dfg) = c("perc_reads", "rootReads", "totalCounts", 
                    "minimizers", "distinctMinimizers", "rank",
                    "taxid", "taxonomy", "sample", "group")
write.table(dfg, paste0(outfile_prefix),
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)