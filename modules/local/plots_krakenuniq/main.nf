process PLOTS_KRAKENUNIQ {

    tag "plots krakenUniq"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table
    val min_reads
    val contaminants

    output:
    path 'krakenuniq_summary_kingdoms.tsv',                  emit: kingdoms
    path 'krakenuniq_removedReadsFromPlots.tsv',             emit: removedReadsFromPlots
    path 'krakenuniq_read_table_at_species_level_postCleaning.tsv', emit: clean_reads
    path '*.tsv'
    path '*.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_krakenuniq.r")

    """
    Rscript $plot_heatmap_script $count_table $min_reads "$contaminants"
    """
}
//  "$contaminants" "$tax_level"
