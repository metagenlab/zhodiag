process KRAKENUNIQ_COMBINE_REPORTS {
    tag "combining krakenUniq reports"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path k2reports
    path(metadata)

    output:
    path "krakenuniq_combined.tsv", emit: combine_long

    script:
    def kraken_merge_script_long = workflow.projectDir.resolve("bin/combine_krakenuniq_reports.r")

    """
    Rscript $kraken_merge_script_long ${metadata} krakenuniq_combined.tsv ${k2reports.join(' ')} 
    """
}
