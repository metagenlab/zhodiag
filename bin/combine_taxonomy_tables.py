#!/usr/bin/env python3
import pandas as pd
import sys
import os
from functools import reduce

files = sys.argv[1:]  # all input files from Nextflow

dfs = []
sample_names = []

for f in files:
    sample_name = os.path.basename(f).replace("_taxonomy.tsv", "")
    sample_names.append(sample_name)
    
    df = pd.read_csv(f, sep="\t")
    df['taxid'] = df['taxid'].astype(str)
    df = df.rename(columns={"nReads": sample_name})
    dfs.append(df)

merge_cols = ['taxid', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

df_merged = reduce(lambda left, right: pd.merge(left, right, on=merge_cols, how='outer'), dfs)

df_merged[sample_names] = df_merged[sample_names].fillna(0).astype(int)

df_merged = df_merged.sort_values(by='taxid')

df_merged.to_csv('combined_taxonomy_counts.tsv', sep='\t', index=False)
