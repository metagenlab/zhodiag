process ANALYSIS_COMBINED {

    tag "Combining tables"
    label 'process_medium'

    container 'docker://lbwang/rocker-genome:latest'
    // container "docker://metagenlab/taxonomizr:1.0"

    input:
    path accession_tables
    path taxid_tables
    path(metadata)

    output:
    path 'combined_table_by_accession.tsv'
    path 'combined_table_by_taxid.tsv'
    path 'heatmap_by_taxid_fillCoverage_labelMappedReads.pdf'
    path 'heatmap_nBasesCovered_by_taxid_genusLevel.pdf'

    script:
    def analysis_script = workflow.projectDir.resolve("bin/analysis_combined.r")

    """
    Rscript $analysis_script ${metadata} accession $accession_tables
    Rscript $analysis_script ${metadata} taxid $taxid_tables
    """
}
