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
    path '*.pdf'
    path '*.tsv'

    script:
    def plot_script = workflow.projectDir.resolve("bin/plots_minimap2.r")

    """
    Rscript $plot_script $paf $mapq $coverage $annotation
    """

}