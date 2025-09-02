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
    def counts_out = kraken2_counts.getSimpleName().replaceFirst(/\.tsv$/, '_full_taxonomy.tsv')
    def minimizers_out = kraken2_minimizers.getSimpleName().replaceFirst(/\.tsv$/, '_full_taxonomy.tsv')

    """
    python3 $kraken_taxonomy_script -i $kraken2_counts
    python3 $kraken_taxonomy_script -i $kraken2_minimizers
    """

}
