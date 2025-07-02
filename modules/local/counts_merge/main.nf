process COUNTS_MERGE {

    tag "counts merge"

    input:
    path taxonomy_files

    output:
    path 'combined_taxonomy_counts.tsv'

    script:
    def script_dir = file("bin")
    def kraken_merge_script = script_dir.resolve("combine_taxonomy_tables.py").toString()
    """
    python3 $kraken_merge_script ${taxonomy_files.join(' ')}
    """
}
