process KRAKEN2_KREPORT2KRONA {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(kreport)

    output:
    tuple val(meta), path("*_2krona.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def filename = kreport.getName()
    def base = filename.replaceFirst(/\.txt$/, '')
    def parts = base.split('_')
    def suffix = parts.size() >= 4 ? parts[1..3].join('_') : 'unknown'
    """
    kreport2krona.py \\
        -r ${kreport} \\
        -o ${prefix}_${suffix}_2krona.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kreport2krona.py: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kreport2krona.py: ${VERSION}
    END_VERSIONS
    """
}
