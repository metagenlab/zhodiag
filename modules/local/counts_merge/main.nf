process COUNTS_MERGE {

    tag "counts merge"
    container "docker://metagenlab/biopython:1.79c"
       
    input:
    path taxonomy_files

    output:
    path 'combined_counts*.tsv'      , emit: counts

    script:
    def script_dir = file("bin")
    def kraken_merge_script = script_dir.resolve("combine_taxonomy_tables.py").toString()
    """
    $kraken_merge_script ${taxonomy_files.join(' ')}
    """
}
