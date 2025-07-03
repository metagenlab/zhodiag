import sys
import csv
from collections import Counter

# Check command-line arguments
if len(sys.argv) != 3:
    print(f"Usage: python3 {sys.argv[0]} <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Count occurrences of the last column
group_counts = Counter()

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        tax_id = line.split()[2]
        group_counts[tax_id] += 1

# Write the results to a TSV file
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["taxID", "nReads"])
    for tax_id, count in group_counts.items():
        writer.writerow([tax_id, count])
