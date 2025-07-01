#!/usr/bin/env python3

# custom script to compute the proportion of k-mers assigned to each taxID.
import sys
from collections import defaultdict

def parse_kraken2_line(line):
    # Kraken2 output fields:
    # 0: classification status (C/U)
    # 1: read ID
    # 2: taxid assigned
    # 3: some counts (e.g. 134|151)
    # 4: list of kmer assignments (taxid:count pairs separated by spaces)
    
    fields = line.strip().split('\t')
    if len(fields) < 5:
        return None
    
    kmer_field = fields[4]
    kmer_parts = kmer_field.split()
    
    taxid_counts = defaultdict(int)
    for part in kmer_parts:
        if ':' not in part:
            continue
        taxid, count = part.split(':')
        try:
            count = int(count)
            taxid_counts[taxid] += count
        except ValueError:
            continue
    
    return taxid_counts

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} kraken2_output_file", file=sys.stderr)
        sys.exit(1)
    
    taxid_totals = defaultdict(int)
    total_kmers = 0
    
    with open(sys.argv[1], 'r') as infile:
        for line in infile:
            if line.startswith('C') or line.startswith('U'):
                counts = parse_kraken2_line(line)
                if counts is None:
                    continue
                for taxid, count in counts.items():
                    taxid_totals[taxid] += count
                    total_kmers += count
    
    # Print report sorted by descending count
    print(f"{'TaxID':<10}\tCount\tPercentage")
    for taxid, count in sorted(taxid_totals.items(), key=lambda x: x[1], reverse=True):
        pct = (count / total_kmers) * 100 if total_kmers > 0 else 0
        print(f"{taxid:<10}\t{count}\t{pct:.2f}%")

if __name__ == "__main__":
    main()
