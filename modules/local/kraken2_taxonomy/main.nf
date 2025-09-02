process KRAKEN2_TAXONOMY {

    tag "Adding full taxonomy"

    conda "${moduleDir}/environment.yml"
    // container './modules/local/kraken2_taxonomy/ete3_custom.sif'
    container "docker://metagenlab/ete:1.0"

    input:
    path kraken2_counts
    path kraken2_minimizers

    output:
    path '*_full_taxonomy.tsv'

    script:
    def script_dir = file("bin")
    def kraken_taxonomy_script = script_dir.resolve("kraken2taxonomy.py").toString()
    def counts_files = kraken2_counts.collect { it.getName() }.join(' ')
    def minimizer_files = kraken2_minimizers.collect { it.getName() }.join(' ')

    """
    for file in $counts_files; do
        python3 $kraken_taxonomy_script -i $kraken2_counts
    done

    for file in $minimizer_files; do
        python3 $kraken_taxonomy_script -i $kraken2_minimizers
    done
    """

}
