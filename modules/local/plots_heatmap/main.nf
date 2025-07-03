process PLOTS_HEATMAP {

    tag "plots heatmap"
    container 'docker://rocker/tidyverse:4.1.0'

    input:
    path count_table

    output:
    path '*heatmap_nTop100*.pdf'
    path '*heatmap_nTop100_exclHumanAmbigUnclas*.pdf'

    script:
    def script_dir = file("bin")
    def plot_heatmap_script = script_dir.resolve("plot_nReads_taxa_heatmap.r").toString()
    def infileName = count_table.getName().replaceFirst(/(\.[^.]+)?$/, '')  // remove extension
    def outfile_prefix = infileName

    """
    Rscript $plot_heatmap_script $count_table $outfile_prefix
    """
}
