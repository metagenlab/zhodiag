#!/usr/bin/env nextflow

process FASTQC {

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
}
