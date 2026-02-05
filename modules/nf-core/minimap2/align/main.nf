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
    val paf_output

    output:
    tuple val(meta), path("*.paf")                  , optional: true, emit: paf
    tuple val(meta), path("*.bam")                  , optional: true, emit: bam
    tuple val(meta), path("*.sam")                  , optional: true, emit: sam
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

    def sam_file = "${prefix}_${genome}.sam"
    def bam_file = "${prefix}_${genome}.bam"
    def paf_file = "${prefix}_${genome}.paf"
    def flagstat_file = "${prefix}_${genome}.flagstat.txt"
    def log_file = "${prefix}_${genome}.log"

    def cigar_paf = cigar_paf_format ? "-c" : ''
    def set_cigar_bam = cigar_bam ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    def read1 = reads[0]

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Determine if read1 is gzipped or not
    if [[ "${read1}" == *.gz ]]; then
      line_count=\$(zcat "$read1" | wc -l)
    else
      line_count=\$(cat "$read1" | wc -l)
    fi

    if [[ \$line_count -eq 0 ]]; then
      echo "Read file is empty. Creating empty output files..."
      touch ${flagstat_file}


      exit 0
    fi

    echo "Running minimap2 alignment..." >&2

    # 1. Align reads to reference, output to SAM
    {
      ln -s ${target} local_genome
      $samtools_reset_fastq \\
      minimap2 -ax sr \\
          --split-prefix ${prefix}_${genome}.sam \\
          $args \\
          -t $task.cpus \\
          $cigar_paf \\
          $set_cigar_bam \\
          local_genome \\
          $query \\
          > ${sam_file}
    } 2> ${log_file}

    ${bam_format ? "samtools view -@ ${task.cpus-1} -b ${sam_file} > ${bam_file}" : ""}

    samtools flagstat ${bam_format ? bam_file : sam_file} > ${flagstat_file}

    ${paf_output ? "paftools.js sam2paf ${sam_file} > ${paf_file}" : ""}

    ${bam_format && unmapped_fq ? """
    samtools view -b -f 1 -F 2 -F 256 ${bam_file} | samtools fastq -1 ${prefix}_${genome}_unmapped_R1.fastq.gz -2 ${prefix}_${genome}_unmapped_R2.fastq.gz /dev/stdin
    """ : ""}


    ${bam_format ? "rm -f ${sam_file}" : ""}

    # Write versions.yml
    minimap2_ver=\$(minimap2 --version 2>&1)
    samtools_ver=\$(samtools --version | head -n1 | sed 's/^.*samtools //')

    cat <<EOF > versions.yml
    minimap2: \$minimap2_ver
    samtools: \$samtools_ver
    EOF
    """
}
