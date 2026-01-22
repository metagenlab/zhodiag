process PLOTS_KRAKENUNIQ {

    tag "plots krakenUniq"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table
    val min_reads
    // val contaminants
    // val tax_level

    output:
    // path 'table_at*totalCounts.tsv',           emit: counts
    // path 'table_at*distinctMinimizers.tsv',    emit: minimizers
    path '*.tsv'
    path '*.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_krakenuniq.r")

    """
    Rscript $plot_heatmap_script $count_table $min_reads
    """
}
//  "$contaminants" "$tax_level"
