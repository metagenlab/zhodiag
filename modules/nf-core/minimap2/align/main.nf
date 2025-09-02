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
    def read2 = reads.size() > 1 ? reads[1] : null
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Check if reads are empty
    if [[ $(zcat ${read1} | wc -l) -eq 0 ]]; then
        echo "Input read file ${reads[0]} is empty. Skipping minimap2." >&2

        # Create empty outputs to avoid Nextflow errors
        ${bam_format ? "touch ${prefix}_${genome}.bam" : ""}
        ${paf_output ? "touch ${prefix}_${genome}.paf" : ""}
        echo -e "0 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    0 + 0 mapped (0.00% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (0.00% : N/A)
    0 + 0 with itself and mate mapped
    0 + 0 singletons (0.00% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)" > ${prefix}_${genome}.flagstat.txt

        ${bam_format && unmapped_fq ? "touch ${prefix}_${genome}_unmapped_R1.fastq.gz ${prefix}_${genome}_unmapped_R2.fastq.gz" : ""}

    else

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

        # Check for @SQ lines in SAM header (i.e., mapped reads)
        if grep -q '^@SQ' ${sam_file}; then
            echo "Mapped reads detected. Continuing with SAM processing..." >&2

            ${bam_format ? "samtools view -@ ${task.cpus-1} -b ${sam_file} > ${bam_file}" : ""}

            samtools flagstat ${bam_format ? bam_file : sam_file} > ${flagstat_file}

            ${paf_output ? "paftools.js sam2paf ${sam_file} > ${paf_file}" : ""}

            ${bam_format && unmapped_fq ? "samtools fastq -f 4 -@ ${task.cpus-1} ${bam_file} -1 ${prefix}_${genome}_unmapped_R1.fastq.gz -2 ${prefix}_${genome}_unmapped_R2.fastq.gz" : ""}

        else
            echo "No mapped reads. Creating empty outputs..." >&2

            ${bam_format ? "samtools view -H ${sam_file} | samtools view -b - > ${bam_file}" : ""}

            echo -e "0 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    0 + 0 mapped (0.00% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (0.00% : N/A)
    0 + 0 with itself and mate mapped
    0 + 0 singletons (0.00% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)" > ${flagstat_file}

            ${paf_output ? "touch ${paf_file}" : ""}
            ${bam_format && unmapped_fq ? "touch ${prefix}_${genome}_unmapped_R1.fastq.gz ${prefix}_${genome}_unmapped_R2.fastq.gz" : ""}
        fi

        ${bam_format ? "rm -f ${sam_file}" : ""}
    fi

    # 7. Write versions
    minimap2_ver=\$(minimap2 --version 2>&1)
    samtools_ver=\$(samtools --version | head -n1 | sed 's/^.*samtools //')

    cat <<EOF > versions.yml
    MINIMAP_HOST:
        minimap2: \$minimap2_ver
        samtools: \$samtools_ver
    EOF
    """


}
