process PLOTS_MAPPING {

    tag "Combining tables"
    label 'process_medium'

    container 'docker://lbwang/rocker-genome:latest'
    // container "docker://metagenlab/taxonomizr:1.0"

    input:
    path(metadata)
    val min_reads
    val contaminants
    path accession_tables
    path taxid_tables

    output:
    path 'bowtie2_summary_kingdoms.tsv',                  emit: kingdoms
    path 'bowtie2_removedReadsFromPlots.tsv',             emit: removedReadsFromPlots
    path 'bowtie2_read_table_at_species_level_postCleaning.tsv', emit: clean_reads
    path '*.pdf'
    path '*.tsv'
    // path 'heatmap_nBasesCovered_by_taxid_genusLevel.pdf'

    script:
    def analysis_script = workflow.projectDir.resolve("bin/plots_mapping.r")

    """
    Rscript $analysis_script ${metadata} accession $min_reads "$contaminants" $accession_tables
    Rscript $analysis_script ${metadata} taxid $min_reads "$contaminants" $taxid_tables
    """
}
