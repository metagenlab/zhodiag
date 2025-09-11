process SUMMARY_STATS_CANDIDATES {
    tag "$meta.id"
    label 'process_medium'

    // container 'docker://lbwang/rocker-genome:latest'
    container "docker://metagenlab/taxonomizr:1.0"

    input:
    tuple val(meta), path(depth_file), path(map_file)

    output:
    tuple val(meta), path("*by_accession.tsv"), emit: accession_table
    tuple val(meta), path("*by_taxid.tsv"), emit: taxid_table
    tuple val(meta), path("*.pdf")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def analysis_script = workflow.projectDir.resolve("bin/analysis_by_sample.r")

    """
    export TAXONOMY_DB=${workflow.projectDir}/data/accessionTaxa.sql
    Rscript $analysis_script $depth_file $map_file ${prefix}
    """
}
