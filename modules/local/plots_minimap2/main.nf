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
    path '*/*/*_boxplot_all.pdf'
    path '*/*/*_heatmap_all.pdf'
    path '*/*/*_filtered_hits.tsv'
    path '*/*/*group_heatmap_all.pdf'
    path '*/*/*group_boxplot_all.pdf'

    script:
    def plot_script = workflow.projectDir.resolve("bin/plots_minimap2.r")
    def out_dir = "${kingdom}/mapq${mapq}"
    def prefix = "${kingdom}_mapq${mapq}"

    """
    mkdir -p ${out_dir}
    Rscript $plot_script $paf ${prefix} $mapq $coverage $annotation

    mv *_heatmap_all.pdf ${out_dir}/
    mv *_boxplot_all.pdf ${out_dir}/
    mv *_filtered_hits.tsv ${out_dir}/
    mv *_group_heatmap_all.pdf ${out_dir}/
    mv *_group_boxplot_all.pdf ${out_dir}/
    """

}