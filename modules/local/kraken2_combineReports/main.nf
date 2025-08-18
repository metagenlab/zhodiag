process KRAKEN2_COMBINE_REPORTS {

    tag "Combining ${kraken_reports.size()} Kraken2 reports"

    input:
    path kraken_reports
    path metadata_table

    output:
    // path "*/*/kraken2_reports_combined_wide*.tsv", emit: combine_wide
    path "kraken2_combined.tsv", emit: combine_long
    path "versions_krakentools.yml", emit: versions


    script:
    // def script_dir = file("bin")
    def kraken_merge_script_long = workflow.projectDir.resolve("bin/combine_kraken2_reports_long.py")

    """
    $kraken_merge_script_long -m ${metadata_table} -i ${kraken_reports} -o kraken2_combined.tsv
    
    
    cat <<-END_VERSIONS > versions_krakentools.yml
    "${task.process}":
        merge_scripts: custom
    END_VERSIONS
"""
}
    // $kraken_merge_script_wide ${reportFiles} ${out_dir}/kraken2_reports_combined_wide_${db}_${conf}.tsv
