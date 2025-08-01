process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    input:
    tuple val(meta), path(reads)
    path reference
    val bam_format
    val cigar_paf_format
    val cigar_bam
    val unmapped_fq

    output:
    tuple val(meta), path("*.paf")                  , optional: true, emit: paf
    tuple val(meta), path("*.bam")                  , optional: true, emit: bam
    tuple val(meta), path("*.flagstat.txt")         , optional: true, emit: flagstat
    tuple val(meta), path("*.fastq.gz")             , optional: true, emit: unmapped
    tuple val(meta), path("*.log")                  , optional: true, emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = reference.getBaseName().replaceFirst(/(\.fna|\.fa|\.fasta)?(\.gz)?$/, '')
    def bam_output_path = "${prefix}_${genome}.bam"
    def bam_output = bam_format ? "-a | samtools view -@ ${task.cpus-1} -b -o ${prefix}_${genome}.bam" : ''
    def logfiles = "${prefix}_${genome}.log"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)
    def flagstat_file = bam_format ? "${prefix}_${genome}.flagstat.txt" : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    {
        $samtools_reset_fastq \\
        minimap2 -x sr \\
            $args \\
            -t $task.cpus \\
            $target \\
            $query \\
            $cigar_paf \\
            $set_cigar_bam \\
            $bam_output
    } 2> $logfiles

    ${bam_format ? "samtools flagstat ${bam_output_path} > ${flagstat_file}" : ""}
    ${bam_format && unmapped_fq ? "samtools fastq -f 4 -@ ${task.cpus-1} ${bam_output_path} -1 ${prefix}_${genome}_unmapped_R1.fastq.gz -2 ${prefix}_${genome}_unmapped_R2.fastq.gz" : ""}

    minimap2_ver=\$(minimap2 --version 2>&1)
    samtools_ver=\$(samtools --version | head -n1 | sed 's/^.*samtools //')

    cat <<EOF > versions.yml
    MINIMAP_HOST:
        minimap2: \$minimap2_ver
        samtools: \$samtools_ver
    EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    touch $output_file
    ${bam_index}

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
  minimap2: \$(minimap2 --version 2>&1)
END_VERSIONS
    """
}
