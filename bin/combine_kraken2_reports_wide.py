#!/usr/bin/env python3
import sys
import os
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
            total_counts = int(cols[2])
            distinct_minimizers = int(cols[4])
            rank_code = cols[5]  # ← correct column for rank
            taxid = cols[6]
            taxonomy = cols[7].strip()
            data[taxid] = (total_counts, distinct_minimizers, rank_code, taxonomy)
    return data

def main(files, output_filename):
    sample_data = defaultdict(dict)  # sample -> taxid -> (totalCounts, distinctMinimizers)
    taxid_to_taxonomy = {}
    taxid_to_rank = {}

    for f in files:
        sample_name = os.path.basename(f).split('_')[0]
        parsed = parse_report(f)
        for taxid, (total_counts, distinct_minimizers, rank_code, taxonomy) in parsed.items():
            sample_data[sample_name][taxid] = (total_counts, distinct_minimizers)
            taxid_to_taxonomy[taxid] = taxonomy
            taxid_to_rank[taxid] = rank_code

    samples = sorted(sample_data.keys())
    all_taxids = set()
    for s in samples:
        all_taxids.update(sample_data[s].keys())

    # Header
    header = ["taxID", "rank", "taxonomy"]
    for s in samples:
        header.append(f"{s}_totalCounts")
        header.append(f"{s}_distinctMinimizers")

    with open(output_filename, "w") as out:
        out.write("\t".join(header) + "\n")

        for taxid in sorted(all_taxids, key=lambda x: int(x) if x.isdigit() else x):
            taxonomy = taxid_to_taxonomy.get(taxid, "")
            rank = taxid_to_rank.get(taxid, "")
            row = [taxid, rank, taxonomy]
            for s in samples:
                counts = sample_data[s].get(taxid, (0, 0))
                row.append(str(counts[0]))
                row.append(str(counts[1]))
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: combine_kraken2_reports.py report1 report2 ... output_filename", file=sys.stderr)
        sys.exit(1)
    report_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    main(report_files, output_file)
