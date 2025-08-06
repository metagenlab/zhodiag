process KRONA_PLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ('*.html'), emit: html
    path "versions_krona.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = report.getName()
    def base = filename.replaceFirst(/\.txt$/, '')
    def parts = base.split('_')
    def suffix = parts.size() >= 4 ? parts[1..2].join('_') : 'unknown'
    """
    ktImportText  \\
        $args \\
        -o ${prefix}_${suffix}.html \\
        $report

    cat <<-END_VERSIONS > versions_krona.yml
    "${task.process}":
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
    END_VERSIONS
    """
}
