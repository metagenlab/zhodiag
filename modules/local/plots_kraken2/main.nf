process PLOTS_KRAKEN2 {

    tag "plots kraken2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table

    output:
    path '*heatmap_all.pdf'
    path '*heatmap_exclHuman.pdf'
    path '*heatmap_exclHuman_facetGroups.pdf'
    path '*totalCounts_vs_distinctMinimizers.pdf'
    path '*heatmap_exclControls.pdf'
    path '*heatmap_exclControls_facetGroups.pdf'
    path '*totalCounts_vs_distinctMinimizers_exclControls.pdf'

    script:
    def script_dir = file("bin")
    def plot_heatmap_script = script_dir.resolve("plots_kraken.r").toString()

    // Extract base name and suffix from file
    def filename = count_table.getName()
    def base = filename.replaceFirst(/\.tsv$/, '')
    def parts = base.split('_')
    def suffix = parts.size() >= 5 ? parts[2..4].join('_') : 'unknown'

    """
    Rscript $plot_heatmap_script $count_table $suffix
    """
}
