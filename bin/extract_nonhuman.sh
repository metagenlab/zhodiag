#!/bin/bash

# Usage: ./extract_nonhuman.sh processed_file.tsv original_kraken2.tsv > filtered_output.tsv

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 processed_file.tsv original_kraken2.tsv"
    exit 1
fi

processed="$1"
original="$2"

# Step 1: Get readIDs to keep
awk -F'\t' '
{
    readid = $1
    decision = $NF

    # Skip if decision is 0 or 9606
    if (decision == "0" || decision == "9606") next

    # If decision is "check", inspect taxIDs
    if (decision == "check") {
        split($2, parts, " ")
        keep = 0
        for (i in parts) {
            split(parts[i], kv, ":")
            taxid = kv[1]
            if (taxid != "0" && taxid != "9606" && taxid != "|") {
                keep = 1
                break
            }
        }
        if (!keep) next
    }

    # Passed filters
    print readid
}
' "$processed" > readIDs_to_keep.txt

# Step 2: Extract original Kraken2 lines where read ID (column 2) matches
awk 'NR==FNR { keep[$1]; next } $2 in keep' readIDs_to_keep.txt "$original"
