process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://quay.io/biocontainers/bbmap:39.26--he5f24ec_0' :
    'docker://quay.io/biocontainers/bbmap:39.26--he5f24ec_0' }"

    input:
    tuple val(meta), path(fastq)
    path ref

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.fastq.gz"), emit: unmapped
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*_bbmap_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def stats = "statsfile=${prefix}_bbmap_stats.txt"
    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"
    def unmapped_output = meta.single_end ?
        "outu=${prefix}_unmapped.fastq.gz" :
        "outu1=${prefix}_unmapped_1.fastq.gz outu2=${prefix}_unmapped_2.fastq.gz"

    // Determine db setting
    def db = ''
    if ( ref.isDirectory() ) {
        if ( ref ==~ /(.\/)?ref\/?/ ) {
            db = ''
        } else {
            db = "path=${ref}"
        }
    } else {
        db = "ref=${ref}"
    }

    """
    set -x
    bbmap.sh \
        $db \
        $input \
        out=stdout.sam \
        $unmapped_output \
        $args \
        $stats \
        threads=$task.cpus \
        -Xmx${task.memory.toGiga()}g | \
        samtools view -Sb - > ${prefix}.bam \
        2> ${prefix}.bbmap.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bbmap.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
