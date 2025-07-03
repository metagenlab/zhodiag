process SORT_INDEX_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/samtools:latest"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${sorted_prefix}.bam"), emit: sorted_bam
    tuple val(meta), path("${sorted_prefix}.bam.bai"), emit: bam_index

    script:
    def sorted_prefix = bam.getBaseName().replaceFirst(/\.bam$/, '') + "_sorted"

    """
    samtools sort -@ ${task.cpus - 1} -o ${sorted_prefix}.bam ${bam}
    samtools index -@ ${task.cpus - 1} ${sorted_prefix}.bam
    """
}