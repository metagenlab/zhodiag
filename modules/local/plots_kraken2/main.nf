process PLOTS_KRAKEN2 {

    tag "plots kraken2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table
    val contaminants
    val tax_level

    output:
    path 'table_at*totalCounts.tsv',           emit: counts
    path 'table_at*distinctMinimizers.tsv',    emit: minimizers
    path '*.tsv'
    path '*.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_kraken.r")

    """
    Rscript $plot_heatmap_script $count_table "$contaminants" "$tax_level"
    """
}
