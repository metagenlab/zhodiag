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
    print(f)
    s = sub("\\.krakenuniq\\.report.*", "", f)
    print(s)
    a = fread(f, header = TRUE, sep = '\t', strip.white = TRUE, fill = TRUE, skip = "taxReads")
    b = a %>% mutate(sample = s)
    df = bind_rows(df, b)
}
df$sample = as.character(df$sample)
dfg = left_join(df, metadata, by = "sample")
dfg = dfg %>% mutate(taxName = gsub("#", "", taxName))

write.table(dfg, paste0(outfile_prefix),
            col.names = TRUE, row.names = FALSE,
            sep = '\t', quote = FALSE)