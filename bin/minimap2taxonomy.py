#!/usr/bin/env python3
"""
Annotate a Kraken2-like TSV with full taxonomy using ete3 NCBITaxa.
Extracts TaxID from target_name column and adds domain, superkingdom, kingdom, phylum, class, order, family, genus, species.
"""

import argparse
import os
import re
import pandas as pd
from ete3 import NCBITaxa

CANONICAL_SUPERKINGDOMS = {
    2759: "Eukaryota",
    2: "Bacteria",
    2157: "Archaea",
    10239: "Viruses",
    12884: "Viroids",
    12908: "Satellites"
}

DEFAULT_RANKS = ["domain","superkingdom","kingdom","phylum","class","order","family","genus","species"]

# Define the expected input columns
INPUT_COLUMNS = ['query_name','query_length','query_start','query_stop','strand',
                 'target_name','target_length','target_start','target_end',
                 'matching_bases','nBases','quality','sample','group']

def get_ncbi_handle():
    """Return dynamic NCBITaxa handle (env override, home dir, fallback)."""
    env_db = os.environ.get("ETE_DBFILE")
    if env_db:
        os.makedirs(os.path.dirname(env_db), exist_ok=True)
        return NCBITaxa(dbfile=env_db)

    home_db = os.path.expanduser("~/.etetoolkit/taxa.sqlite")
    try:
        os.makedirs(os.path.dirname(home_db), exist_ok=True)
        return NCBITaxa(dbfile=home_db)
    except Exception:
        fallback_db = os.path.join(os.getcwd(), ".etetoolkit", "taxa.sqlite")
        os.makedirs(os.path.dirname(fallback_db), exist_ok=True)
        return NCBITaxa(dbfile=fallback_db)

def build_rank_map(ncbi: NCBITaxa, taxid: int, ranks=DEFAULT_RANKS):
    """Return {rank -> name} for a single taxid, filling domain/superkingdom robustly."""
    out = {r: None for r in ranks}
    if taxid <= 1:
        return out
    try:
        lineage = ncbi.get_lineage(taxid)
        if not lineage:
            return out
        lineage_ranks = ncbi.get_rank(lineage)
        lineage_names = ncbi.get_taxid_translator(lineage)
        for tid, r in lineage_ranks.items():
            if r in out:
                out[r] = lineage_names.get(tid)
        this_rank = ncbi.get_rank([taxid]).get(taxid)
        if this_rank in out:
            out[this_rank] = ncbi.get_taxid_translator([taxid])[taxid]
        if out.get("superkingdom") in (None, "", "no rank"):
            for sk_tid, sk_name in CANONICAL_SUPERKINGDOMS.items():
                if sk_tid in lineage:
                    out["superkingdom"] = sk_name
                    break
        if out.get("domain") in (None, "", "no rank"):
            out["domain"] = out.get("superkingdom")
        return out
    except Exception:
        return out

def extract_taxid(target_name: str) -> int:
    """Extract numeric taxid from Kraken-style string."""
    m = re.search(r"kraken:taxid\|(\d+)", target_name)
    if m:
        return int(m.group(1))
    return 0

def main():
    parser = argparse.ArgumentParser(description="Annotate TSV with taxonomy from taxid")
    parser.add_argument("-i", "--input", required=True, help="Input TSV")
    parser.add_argument("-o", "--output", required=True, help="Output TSV with taxonomy columns")
    args = parser.parse_args()

    # Read input with explicit column names
    df = pd.read_csv(args.input, sep="\t", names=INPUT_COLUMNS, dtype=str, header=None)

    # Extract taxid from target_name
    df["taxid"] = df["target_name"].apply(extract_taxid)
    taxids = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")
    unique_taxids = pd.unique(taxids.dropna())

    # NCBI handle
    ncbi = get_ncbi_handle()

    # Build taxonomy maps once per unique taxid
    maps = {tid: build_rank_map(ncbi, tid, ranks=DEFAULT_RANKS) for tid in unique_taxids}
    tax_maps_df = pd.DataFrame.from_dict(maps, orient="index").rename_axis("taxid").reset_index()
    tax_maps_df["taxid"] = pd.to_numeric(tax_maps_df["taxid"], errors="coerce").astype("Int64")

    # Merge taxonomy back to main dataframe
    df_out = df.merge(tax_maps_df, on="taxid", how="left")

    df_out.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
