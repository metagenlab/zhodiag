process KRAKEN2_COMBINE_REPORTS {
    tag "combining kraken2 reports"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path k2reports
    path(metadata)

    output:
    path "kraken2_combined.tsv", emit: combine_long

    script:
    def kraken_merge_script_long = workflow.projectDir.resolve("bin/combine_kraken_reports.r")

    """
    Rscript $kraken_merge_script_long ${metadata} kraken2_combined.tsv ${k2reports.join(' ')} 
    """
}
