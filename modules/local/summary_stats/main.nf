process SUMMARY_STATS_CANDIDATES {

    tag "$meta.id"
    label 'process_medium'

    container 'docker://lbwang/rocker-genome:latest'

    input:
    tuple val(meta), path(input_table)
    tuple val(meta), path(map_summary_table)

    output:
    path '*.tsv'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def analysis_script = workflow.projectDir.resolve("bin/depth_tables_by_sample.r")
    def output_prefix = "${prefix}_"

    """
    Rscript $analysis_script $input_table $map_summary_table $output_prefix
    """
}
