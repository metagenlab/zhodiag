process CUSTOM_STAT_REPORT {

    tag "Statistics Report"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    // path(samples)
    path(multiqc_data)
    val trim_tool
    val mapper
    path(host_fasta)
    val krakenuniq
    path krakenuniq_report
    val kraken2
    path kraken2_report
    val map_classified

    output:
    path '*.tsv'
    path '*.pdf'

    script:
    def custom_report_script = workflow.projectDir.resolve("bin/custom_stat_report.r")
    def host_name = host_fasta.getBaseName().replaceFirst(/(\.fna|\.fa|\.fasta)?(\.gz)?$/, '')

    """
    Rscript $custom_report_script \
            $multiqc_data \
            $trim_tool \
            $mapper \
            $host_name \
            $krakenuniq \
            ${krakenuniq ? krakenuniq_report : ''} \
            $kraken2 \
            ${kraken2 ? kraken2_report : ''} \
            $map_classified
    """

}

