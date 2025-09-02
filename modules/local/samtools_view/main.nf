process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    val quality
    val cov

    output:
    tuple val(meta), path("${prefix}_filtered.bam"),   emit: bam                                 emit: bam,              optional: true
    path  "versions_samtools_view.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def output_file = "${prefix}_filtered.bam"

    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        $args \\
        -q $quality \\
        -e '(sclen/qlen) > $cov' \\
        -o ${output_file} \\
        $input

    cat <<-END_VERSIONS > versions_samtools_view.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
