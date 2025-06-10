#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.28--pyhdfd78af_0' :
        'biocontainers/multiqc:1.28--pyhdfd78af_0' }"

    input:
    path '*' 
    val output_name


    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc . --force -n ${output_name}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

}