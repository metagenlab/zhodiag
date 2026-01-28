process CUSTOM_STAT_REPORT {

    tag "Custom Run Statistics Report"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    // path(samples)
    path(multiqc_data)
    val trim_tool
    val mapper
    val krakenuniq
    path krakenuniq_report
    val kraken2
    path kraken2_report

    output:
    path '*.tsv'
    path '*.pdf'

    script:
    def custom_report_script = workflow.projectDir.resolve("bin/custom_stat_report.r")

    """
    Rscript $custom_report_script \
            $multiqc_data \
            $trim_tool \
            $mapper \
            $krakenuniq \
            ${krakenuniq ? krakenuniq_report : ''} \
            $kraken2 \
            ${kraken2 ? kraken2_report : ''}
    """
}
