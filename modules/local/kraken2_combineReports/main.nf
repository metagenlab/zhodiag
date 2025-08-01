process KRAKEN2_COMBINE_REPORTS {

    tag "Combining ${samples_list.size()} Kraken2 reports"

    input:
    val samples_list

    output:
    path "kraken2_reports_combined_wide.tsv", emit: combine_wide
    path "kraken2_combined*.tsv", emit: combine_long

    script:
    // def script_dir = file("bin")
    def kraken_merge_script_long = workflow.projectDir.resolve("bin/combine_kraken2_reports_long.py")
    def kraken_merge_script_wide = workflow.projectDir.resolve("bin/combine_kraken2_reports_wide.py")


    // Validate and sanitize input list of [id, group, path]
    def sanitized = samples_list.collect { sample ->
        if (sample instanceof List || sample instanceof Tuple) {
            def id = sample[0].toString()
            def group = sample[1].toString()
            def reportPath = sample[2].toString()
            // Check if file exists, fail early if not
            if (!file(reportPath).exists()) {
                throw new RuntimeException("Kraken2 report file not found: ${reportPath}")
            }
            return [id, group, reportPath]
        } else {
            throw new RuntimeException("Invalid sample format: ${sample.dump()}")
        }
    }

    // Prepare metadata lines: sample_id, group_label, report_filename
    def metadata_lines = sanitized.collect { entry ->
        def sid = entry[0]
        def group = entry[1]
        def reportFileName = new File(entry[2]).getName()
        return "${sid}\t${group}\t${reportFileName}"
    }.join('\n')

    def metadataFile = "metadata.txt"
    def reportFiles = sanitized.collect { it[2] }.join(' ')

    """
    echo -e '${metadata_lines}' > ${metadataFile}

    $kraken_merge_script_wide ${reportFiles} kraken2_reports_combined_wide.tsv
    $kraken_merge_script_long ${reportFiles} ${metadataFile} kraken2_combined
    """
}
