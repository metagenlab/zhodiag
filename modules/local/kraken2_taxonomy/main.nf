process KRAKEN2_TAXONOMY {

    tag "Adding full taxonomy"

    conda "${moduleDir}/environment.yml"
    container "docker://metagenlab/ete:1.0"

    input:
    path kraken2_counts
    path kraken2_minimizers

    output:
    path '*_full_taxonomy.tsv'

    script:
    def kraken_taxonomy_script = workflow.projectDir.resolve("kraken2taxonomy.py")

    def count_cmds = kraken2_counts.collect { f -> "python3 $kraken_taxonomy_script -i ${f.getName()}" }.join('\n')
    def minimizer_cmds = kraken2_minimizers.collect { f -> "python3 $kraken_taxonomy_script -i ${f.getName()}" }.join('\n')

    """
    $count_cmds

    $minimizer_cmds
    """
}
