#!/usr/bin/env python3
"""
Annotate a Kraken2-like TSV with full taxonomy using ete3 NCBITaxa.
Extracts TaxID from 'taxid' column and adds domain, superkingdom, kingdom, phylum, class, order, family, genus, species.
"""

import argparse
import os
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

DEFAULT_RANKS = ["domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

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

def get_full_lineage(taxid, ncbi):
    """Return dictionary of taxonomy ranks for a given taxid."""
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        lineage_ranks = {r: names[t] for t, r in ranks.items() if r in DEFAULT_RANKS}
        # Add superkingdom/domain if missing from standard ranks
        if "superkingdom" not in lineage_ranks:
            for tid in lineage:
                if tid in CANONICAL_SUPERKINGDOMS:
                    lineage_ranks["superkingdom"] = CANONICAL_SUPERKINGDOMS[tid]
                    break
        if "domain" not in lineage_ranks and "superkingdom" in lineage_ranks:
            lineage_ranks["domain"] = lineage_ranks["superkingdom"]

        return {rank: lineage_ranks.get(rank, "NA") for rank in DEFAULT_RANKS}
    except Exception as e:
        return {rank: "NA" for rank in DEFAULT_RANKS}

def annotate_table(input_file, output_file=None):
    """Load TSV, annotate with taxonomy, write to new TSV."""
    print(f"📥 Reading input: {input_file}")
    df = pd.read_csv(input_file, sep="\t")

    if 'taxid' not in df.columns:
        raise ValueError("Input file must have a 'taxid' column.")

    ncbi = get_ncbi_handle()

    print("🔍 Fetching taxonomy for all taxids...")
    tax_data = df['taxid'].apply(lambda tid: get_full_lineage(int(tid), ncbi))
    tax_df = pd.DataFrame(tax_data.tolist())

    print("🧬 Merging taxonomy with input table...")
    final_df = pd.concat([df, tax_df[DEFAULT_RANKS]], axis=1)

    if not output_file:
        base, ext = os.path.splitext(input_file)
        output_file = base + "_full_taxonomy.tsv"

    print(f"💾 Writing annotated file to: {output_file}")
    final_df.to_csv(output_file, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Add full taxonomy using ETE3 NCBITaxa to Kraken2-like TSV.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with 'taxid' column.")
    parser.add_argument("-o", "--output", help="Output filename. If not set, adds _full_taxonomy.tsv suffix.")
    args = parser.parse_args()

    annotate_table(args.input, args.output)

if __name__ == "__main__":
    main()
