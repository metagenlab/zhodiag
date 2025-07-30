#!/usr/bin/env python3
import sys
import os
import csv
from collections import defaultdict

def parse_report(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 8:
                continue
            root_reads = int(cols[1])
            total_counts = int(cols[2])
            distinct_minimizers = int(cols[4])
            rank = cols[5]
            taxid = cols[6]
            taxonomy = cols[7].strip()
            data[taxid] = (total_counts, distinct_minimizers, rank, taxonomy, root_reads)
    return data

def read_metadata(metadata_path):
    mapping = {}  # filename → (sample_id, group)
    with open(metadata_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue
            sample_id, group, filename = parts
            mapping[filename] = (sample_id, group)
    return mapping

def main(report_files, metadata_path, output_filename_base):
    # Extract suffix from first filename
    first_file = os.path.basename(report_files[0])
    
    # Remove file extension first (e.g., .report.txt)
    base = first_file.replace('.report.txt', '').replace('.report', '')

    parts = base.split('_')
    if len(parts) < 4:
        suffix = "unknown"
    else:
        suffix = '_'.join(parts[1:4])  # e.g., pluspf_20250402.conf_0.5

    output_filename = f"{output_filename_base.rstrip('.tsv')}_{suffix}.tsv"
        
    metadata = read_metadata(metadata_path)

    rows = []
    seen = set()

    for f in report_files:
        filename = os.path.basename(f)
        if filename not in metadata:
            print(f"Warning: {filename} not found in metadata", file=sys.stderr)
            continue
        sample_id, group = metadata[filename]

        parsed = parse_report(f)
        for taxid, (total_counts, distinct_minimizers, rank, taxonomy, root_reads) in parsed.items():
            key = (taxid, sample_id)
            if key in seen:
                continue
            seen.add(key)
            rows.append([
                taxid,
                taxonomy,
                rank,
                sample_id,
                group,
                total_counts,
                distinct_minimizers,
                root_reads
            ])

    # Sort output
    rows.sort(key=lambda x: (int(x[0]) if x[0].isdigit() else x[0], x[3]))

    with open(output_filename, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "taxID", "taxonomy", "rank", "sample", "group",
            "totalCounts", "distinctMinimizers", "rootReads"
        ])
        writer.writerows(rows)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: combine_kraken2_reports_long.py report1 report2 ... metadata.txt output_file_prefix", file=sys.stderr)
        sys.exit(1)
    report_files = sys.argv[1:-2]
    metadata_file = sys.argv[-2]
    output_file_prefix = sys.argv[-1]
    main(report_files, metadata_file, output_file_prefix)
