process SLIM_SAM2BAM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = "${prefix}.bam"

    """
    # 1. Stream slim SAM: only @SQ lines for mapped contigs + other headers + all alignments
    (
      # Contigs with mapped reads
      samtools view -@ ${task.cpus-1} $input | cut -f3 | sort -u | \
      xargs -I{} grep -E "^@SQ.*SN:{}\$" <(samtools view -H $input)

      # Other header lines
      samtools view -H $input | grep -v '^@SQ'

      # All alignments
      samtools view -@ ${task.cpus-1} $input
    ) > ${prefix}_mapped_slim.sam

    # 2. Convert slim SAM to BAM
    samtools view -@ ${task.cpus-1} -b ${prefix}_mapped_slim.sam > ${bam_file}

    # 3. Clean up slim SAM
    rm -f ${prefix}_mapped_slim.sam
    """
}
