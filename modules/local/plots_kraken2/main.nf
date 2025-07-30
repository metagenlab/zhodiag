process PLOTS_KRAKEN2 {

    tag "plots kraken2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table

    output:
    path '*summary_kingdoms.tsv'
    path '*heatmap_all.pdf'
    // path '*heatmap_exclHuman.pdf'
    // path '*heatmap_exclHuman_facetGroups.pdf'
    // path '*totalCounts_vs_distinctMinimizers.pdf'
    path '*heatmap_filtContam.pdf'
    path '*heatmap_filtContam_facetGroups.pdf'
    path '*totalCounts_vs_distinctMinimizers_filtContam.pdf'
    // path '*boxplot_totalCounts_exclHuman.pdf'
    // path '*boxplot_distinctMinimizers_exclHuman.pdf'
    path '*boxplot_totalCounts_filtContam.pdf'
    path '*boxplot_distinctMinimizers_filtContam.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_kraken.r")

    // Extract base name and suffix from file
    def filename = count_table.getName()
    def base = filename.replaceFirst(/\.tsv$/, '')
    def parts = base.split('_')
    def suffix = parts.size() >= 5 ? parts[2..4].join('_') : 'unknown'

    """
    Rscript $plot_heatmap_script $count_table $suffix
    """
}
