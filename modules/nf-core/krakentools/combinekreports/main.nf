process KRAKEN2_COMBINEKREPORTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    path(kreports)

    output:
    path("krakentools_combine_kreports.tsv"), emit: txt
    path "versions_krakentools.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    combine_kreports.py \\
        -r ${kreports} \\
        -o krakentools_combine_kreports.tsv \\
        --display-headers \\
        ${args}

    cat <<-END_VERSIONS > versions_krakentools.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions_krakentools.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """
}
