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
    # 1. Extract contigs with at least one alignment
    samtools view -@ ${task.cpus-1} $input | cut -f3 | sort -u > ${prefix}_mapped_refs.txt

    # 2. Build slim header (only mapped contigs)
    samtools view -H $input | grep -E '^@SQ' | grep -Ff ${prefix}_mapped_refs.txt > ${prefix}_header.sam
    samtools view -H $input | grep -v '^@SQ' >> ${prefix}_header.sam

    # 3. Concatenate slim header + all alignments
    samtools view -@ ${task.cpus-1} $input > ${prefix}_body.sam
    cat ${prefix}_header.sam ${prefix}_body.sam > ${prefix}_mapped_slim.sam

    # 4. Convert to BAM
    samtools view -@ ${task.cpus-1} -b ${prefix}_mapped_slim.sam > ${output_file}

    # 5. Clean up intermediates
    rm -f ${prefix}_header.sam ${prefix}_body.sam ${prefix}_mapped_slim.sam ${prefix}_mapped_refs.txt
    """
}
