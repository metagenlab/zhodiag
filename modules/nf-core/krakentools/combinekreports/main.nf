process KRAKEN2_COMBINEKREPORTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    path(kreports)

    output:
    path("*/*/krakentools_combine_kreports.tsv"), emit: txt
    path("*/*/versions_krakentools.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.2'

    def first_report = kreports instanceof List ? kreports[0] : kreports
    def filename = first_report.getName()  // e.g., 101999248520_pluspf20250402_conf0.5_report.txt

    def matcher = filename =~ /^.+?_(.+?)_(conf\d+(?:\.\d+)?)_report\.txt$/
    if (!matcher.matches()) {
        throw new RuntimeException("Filename does not match expected format: ${filename}")
    }

    def db   = matcher[0][1]  // pluspf20250402
    def conf = matcher[0][2]  // conf0.5
    def out_dir = "${db}/${conf}"

    def output_file = "${out_dir}/krakentools_combine_kreports.tsv"
    def version_file = "${out_dir}/versions_krakentools.yml"

    """
    mkdir -p ${out_dir}

    combine_kreports.py \\
        -r ${kreports.collect { it.getName() }.join(' ')} \\
        -o ${output_file} \\
        --display-headers \\
        ${args}

    cat <<-END_VERSIONS > ${version_file}
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """
}
