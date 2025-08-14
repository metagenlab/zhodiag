process PLOTS_KRAKEN2 {

    tag "plots kraken2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path count_table
    val contaminants

    output:
    path '*/*/*summary_kingdoms.tsv'
    path '*/*/*all_heatmap_totalCounts.pdf'
    path '*/*/*all_heatmap_distinctMinimizers.pdf'
    path '*/*/*group_heatmap_totalCounts.pdf'
    // path '*heatmap_exclHuman.pdf'
    // path '*heatmap_exclHuman_facetGroups.pdf'
    // path '*totalCounts_vs_distinctMinimizers.pdf'
    // path '*heatmap_filtContam.pdf'
    // path '*/*/*heatmap_totalCounts_filtContam_facetGroups.pdf'
    // path '*/*/*heatmap_distinctMinimizers_filtContam_facetGroups.pdf'
    // path '*/*/*group_heatmap_totalCounts_filtContam.pdf'
    path '*/*/*totalCounts_vs_distinctMinimizers.pdf'
    // path '*boxplot_totalCounts_exclHuman.pdf'
    // path '*boxplot_distinctMinimizers_exclHuman.pdf'
    // path '*/*/*filtContam_boxplot_totalCounts.pdf'
    // path '*/*/*filtContam_boxplot_distinctMinimizers.pdf'
    // path '*/*/*group_boxplot_totalCounts_filtContam.pdf'

    script:
    def plot_heatmap_script = workflow.projectDir.resolve("bin/plots_kraken.r")

    def count_table_name = count_table.getName()
    def filename = count_table.getName()

    def matcher = filename =~ /^kraken2_combined_(.+?)_(conf\d+(?:\.\d+)?)_report\.tsv$/
    if (!matcher.matches()) {
        throw new RuntimeException("Filename does not match expected pattern: ${filename}")
    }

    def db = matcher[0][1]
    def conf = matcher[0][2]
    def out_dir = "${db}/${conf}"

    // // Extract base name and suffix from file
    // def filename = count_table.getName()
    // def base = filename.replaceFirst(/\.tsv$/, '')
    // def parts = base.split('_')
    // def suffix = parts.size() >= 5 ? parts[2..4].join('_') : 'unknown'

    """
    mkdir -p ${out_dir}
    Rscript $plot_heatmap_script $count_table ${db}_${conf} $contaminants

    mv *summary_kingdoms.tsv ${out_dir}/
    mv *all_heatmap_totalCounts.pdf ${out_dir}/
    mv *all_heatmap_distinctMinimizers.pdf ${out_dir}/
    mv *group_heatmap_totalCounts.pdf ${out_dir}/
    mv *totalCounts_vs_distinctMinimizers.pdf ${out_dir}/
    """
}
