process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    // tuple val(meta2), path(intervals)

    output:
    // tuple val(meta), path("*depth.tsv"), emit: depth
    tuple val(meta), path("*nonHuman.tsv"), emit: depth_nonHuman
    path "versions_samtoolsdepth.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output =  "${prefix}_depth.tsv"
    def output_nonHuman = "${prefix}_depth_nonHuman.tsv"
    // def positions = intervals ? "-b ${intervals}" : ""
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        -o ${output} \\
        $input
    
    awk -F'|' '\$2 !~ /^9606([^0-9]|\$)/' "${output}" > "${output_nonHuman}"

    cat <<-END_VERSIONS > versions_samtoolsdepth.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
