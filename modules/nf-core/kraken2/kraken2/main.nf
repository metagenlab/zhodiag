process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_veryhigh'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/29ed8f68315625eca61a3de9fcb7b8739fe8da23f5779eda3792b9d276aa3b8f/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:45764814c4bb5bf3' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val confidence
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path("*_classified{.,_}*"), optional: true, emit: classified_reads_fastq
    tuple val(meta), path("*_unclassified{.,_}*"), optional: true, emit: unclassified_reads_fastq
    tuple val(meta), path("*_classifiedreads.txt"), optional: true, emit: classified_reads_assignment
    tuple val(meta), path("*_report.txt"), emit: report
    path "versions_kraken2.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def paired = meta.single_end ? "" : "--paired"
    def conf = confidence
    def out_dir = "."
    def prefix = task.ext.prefix ?: "${meta.id}"
    def classified   = meta.single_end ? "${prefix}_classified.fastq"   : "${prefix}_classified_#.fastq"
    def unclassified = meta.single_end ? "${prefix}_unclassified.fastq" : "${prefix}_unclassified_#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}_classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -f -p $task.cpus *.fastq" : ""

    """
    mkdir -p ${out_dir}

    kraken2 \\
        --db $db \\
        --confidence $conf \\
        --threads $task.cpus \\
        --report ${prefix}_report.txt \\
        --report-minimizer-data \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        $reads


    ${compress_reads_command}

    cat <<-END_VERSIONS > versions_kraken2.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS

    """
}