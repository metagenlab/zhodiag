process COLLECT_REPORTS {
    tag "multiqc_input"

    input:
    path(report_files)

    output:
    path("multiqc_input"), emit: dir

    script:
    """
    mkdir multiqc_input
    cp -r ${report_files.join(' ')} multiqc_input/
    """
}
