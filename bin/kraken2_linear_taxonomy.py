#!/usr/bin/env python

import pandas as pd
from ete3 import NCBITaxa
import argparse

def get_lineage(taxid, main_ranks, ncbi):
    row = {rank: '' for rank in main_ranks}
    
    # Handle known special cases first
    if str(taxid) == "0":
        row = {rank: 'unclassified' for rank in main_ranks}
        del row["domain"]
        row['taxid'] = 0
        return row
    if taxid == "ambiguous":
        row = {rank: 'ambiguous' for rank in main_ranks}
        del row["domain"]
        row['taxid'] = "ambiguous"
        return row
    
    # Now check if taxid is a valid integer
    try:
        taxid_int = int(taxid)
    except ValueError:
        # taxid is not a number, so return NAs with original taxid string
        row = {rank: 'NA' for rank in main_ranks}
        del row["domain"]
        row['taxid'] = taxid
        return row
    
    # taxid is valid integer, do normal lineage lookup
    lineage = ncbi.get_lineage(taxid_int)
    rank_dict = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)

    row['taxid'] = taxid_int
    for lineage_taxid in lineage:
        rank = rank_dict.get(lineage_taxid, '')
        if rank in main_ranks:
            row[rank] = names.get(lineage_taxid, '')

    return row

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="input taxid to counts")
    parser.add_argument("-o", "--output", type=str, help="output", default="taxonomy.tsv")
    parser.add_argument("--db", required=True, help="Path to ete3 taxa.sqlite database")
    args = parser.parse_args()

    # Initialize NCBITaxa *after* parsing arguments
    ncbi = NCBITaxa(dbfile=args.db)

    df = pd.read_csv(args.input, sep="\t")
    taxids = df["taxID"].unique()

    main_ranks = ['acellular root', 'domain', 'kingdom',
                  'phylum', 'class', 'order', 'family', 'genus', 'species']

    rows = [get_lineage(i, main_ranks, ncbi) for i in taxids]
    df_taxo = pd.DataFrame(rows).set_index("taxid")

    df_taxo["superkingdom"] = df_taxo["acellular root"].fillna("") + df_taxo["domain"].fillna("")
    del df_taxo["acellular root"]
    del df_taxo["domain"]

    ordered_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # Convert indices to string for consistent merging
    df_taxo.index = df_taxo.index.astype(str)
    df_counts = df.set_index("taxID")
    df_counts.index = df_counts.index.astype(str)

    # Merge taxonomy with counts (nReads)
    df_out = df_taxo.join(df_counts[["nReads"]], how="left")

    # Save merged output (taxonomy with nReads)
    df_out[ordered_ranks + ["nReads"]].to_csv(args.output, sep="\t")
