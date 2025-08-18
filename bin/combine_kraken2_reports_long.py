#!/usr/bin/env python3
import sys
import os
import csv
import pandas 
import re


def main(report_files, metadata_path, output_filename):
    # Extract suffix from first filename
    first_file = os.path.basename(report_files[0])
    
    # Remove file extension first (e.g., .report.txt)
    base = first_file.replace('.report.txt', '').replace('.report', '')

    parts = base.split('_')
    if len(parts) < 4:
        suffix = "unknown"
    else:
        suffix = '_'.join(parts[1:4])  # e.g., pluspf_20250402.conf_0.5


    df_metadata = pandas.read_csv(metadata_path)[["group"]]
  
    all_df = []

    for f in report_files:

        filename = os.path.basename(f)
        # 1019666478_pandb_conf0.55_report.txt
        m = re.match("(.+)_(.+)_conf([\.\d]+)_report.txt", filename)

        sample, db, conf = m.groups()
        conf = f"conf{conf}"

        df = pandas.read_csv(f, sep="\t", names=["percent", 
                                                 "rootReads", 
                                                 "totalCounts", 
                                                 "totalMinimizers", 
                                                 "distinctMinimizers", 
                                                 "rank", 
                                                 "taxID", 
                                                 "taxonomy"])
        df["sample"] = sample
        all_df.append(df)

    df_merged = pandas.concat(all_df)

    df_merged = df_merged.set_index("sample").join(df_metadata).reset_index()
    df_merged = df_merged.rename(columns={"index": "sample"})

    df_merged.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", '--report_files',
        nargs="+",
        help="Kraken report file.",
    )

    parser.add_argument(
        "-m", '--metadata',
        type=Path,
        help="Samtools depth file (TSV format).",
    )
    parser.add_argument(
        "-o", '--output_prefix',
        type=str,
        default='out',
        help="Output prefix",
    )

    args = parser.parse_args()

    main(args.report_files, args.metadata, args.output_prefix)
