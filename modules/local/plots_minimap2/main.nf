process PLOTS_MINIMAP2 {
    tag "plots minimap2"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path(paf)
    val mapq
    val coverage
    val kingdom
    path(annotation)

    output:
    path '*/*/*_all_boxplot.pdf'
    path '*/*/*_all_heatmap.pdf'
    path '*/*/*_filtered_hits.tsv'
    path '*/*/*group_heatmap.pdf'
    path '*/*/*group_boxplot.pdf'

    script:
    def plot_script = workflow.projectDir.resolve("bin/plots_minimap2.r")
    def out_dir = "${kingdom}/mapq${mapq}"
    def prefix = "${kingdom}_mapq${mapq}"

    """
    mkdir -p ${out_dir}
    Rscript $plot_script $paf ${prefix} $mapq $coverage $annotation

    mv *_all_heatmap.pdf ${out_dir}/
    mv *_all_boxplot.pdf ${out_dir}/
    mv *_filtered_hits.tsv ${out_dir}/
    mv *_group_heatmap.pdf ${out_dir}/
    mv *_group_boxplot.pdf ${out_dir}/
    """

}