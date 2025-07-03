#!/opt/conda/bin/python
import pandas as pd
import sys
import os
from functools import reduce

files = sys.argv[1:]  # all input files from Nextflow

dfs = []
sample_names = []

for f in files:
    sample_name = os.path.basename(f).split('_')[0]
    sample_names.append(sample_name)
    
    df = pd.read_csv(f, sep="\t")
    df['taxid'] = df['taxid'].astype(str)
    df = df.rename(columns={"nReads": sample_name})
    dfs.append(df)

merge_cols = ['taxid', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

df_merged = reduce(lambda left, right: pd.merge(left, right, on=merge_cols, how='outer'), dfs)

df_merged[sample_names] = df_merged[sample_names].fillna(0).astype(int)

df_merged = df_merged.sort_values(by='taxid')

taxonomy_levels = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
suffixes = {
    'kingdom': ' k.',
    'phylum': ' p.',
    'class': ' c.',
    'order': ' o.',
    'family': ' f.',
    'genus': ' g.',
    'species': ' s.'
}

def recursive_fill(row):
    last_valid = None
    for level in taxonomy_levels:
        val = row[level]
        if pd.isna(val) or val == "":
            if last_valid:
                suffix = suffixes.get(level, '')
                row[level] = last_valid + suffix
        else:
            last_valid = val
    return row

df_merged = df_merged.apply(recursive_fill, axis=1)

# Extract suffix from the first filename
first_file = os.path.basename(files[0])
parts = first_file.split('_')
# Get last two parts joined with underscore (for example: "postReassign_0.65.tsv")
suffix_part = "_" + parts[-3] + "_" + parts[-2].replace(".tsv", "")

output_filename = f"combined_counts{suffix_part}.tsv"
df_merged.to_csv(output_filename, sep='\t', index=False)