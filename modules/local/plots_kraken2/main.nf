process PLOTS_KRAKEN2 {

    tag "plots kraken2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table
    val min_reads
    val contaminants

    output:
    path 'kraken2_summary_kingdoms.tsv',               emit: kingdoms
    path 'kraken2_removedReadsFromPlots.tsv',          emit: removedReadsFromPlots
    path 'kraken2_read_table_at_species_level_postCleaning.tsv', emit: clean_reads
    path '*.tsv'
    path '*.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_kraken.r")

    """
    Rscript $plot_heatmap_script $count_table $min_reads "$contaminants"
    """
}
