#!/bin/bash

# Usage: ./process_kmers.sh input_file.tsv [threshold]
# Default threshold is 0.65

if [ -z "$1" ]; then
    echo "Usage: $0 input_file.tsv [threshold]"
    exit 1
fi

input_file="$1"
threshold="${2:-0.65}"  # Default to 0.65 if not provided

gawk -F'\t' -v threshold="$threshold" '
function sort_descending(arr, keys,   i, j, tmp, n) {
    n = 0
    for (i in arr) {
        keys[++n] = i
    }
    # Bubble sort by value (descending)
    for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (arr[keys[j]] > arr[keys[i]]) {
                tmp = keys[i]
                keys[i] = keys[j]
                keys[j] = tmp
            }
        }
    }
    return n
}

{
    split($5, parts, " ");
    delete ref_counts;
    total = 0;

    for (i in parts) {
        split(parts[i], kv, ":");
        ref = kv[1];
        count = kv[2];
        ref_counts[ref] += count;
        total += count;
    }

    delete frac_map;
    for (ref in ref_counts) {
        frac_map[ref] = ref_counts[ref] / total;
    }

    delete sorted_refs;
    n = sort_descending(frac_map, sorted_refs);

    output = "";
    for (i = 1; i <= n; i++) {
        ref = sorted_refs[i];
        output = output sprintf("%s:%.4f ", ref, frac_map[ref]);
    }

    sub(/[ \t]+$/, "", output);  # Trim trailing space

    top_ref = sorted_refs[1];
    top_score = frac_map[top_ref];

    if (top_score > threshold)
        decision = top_ref;
    else
        decision = "check";

    print $1 "\t" output "\t" decision;
}
' "$input_file"
