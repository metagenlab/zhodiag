process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    val quality
    val cov

    output:
    tuple val(meta), path("*_filtered_sorted.bam"),  emit: filtered
    tuple val(meta), path("*_filtered_sorted.bam.bai"),  emit: index
    tuple val(meta), path("*.flagstat.txt")         , optional: true, emit: flagstat

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def output_file = "${prefix}_filtered_sorted.bam"
    def output_flagstat = "${prefix}_filtered_noHuman.flagstat.tsv"

    """
    # filter SAM by quality (= 0 & >mapq) and remove human hits. Output alignments only in sam
    samtools view --threads ${task.cpus-1} $args $input \\
    | awk 'BEGIN{OFS="\\t"} 
            !/^@/ { 
                if (\$5 > $quality && length(\$10) >= 120) {
                    split(\$3, a, "|"); 
                    if (a[2] != 9606) print 
                } 
            }' \\
    > ${prefix}_alignments.sam

    # extract hit IDs
    cut -f3 ${prefix}_alignments.sam | sort -u > ${prefix}_mapped_refs.txt

    # Build slim header (only mapped identifiers)
    awk 'NR==FNR { keep[\$1]; next }
     /^@SQ/ { if (substr(\$2,4) in keep) print; next }
     { print }' ${prefix}_mapped_refs.txt <(samtools view -H $input) > ${prefix}_header.sam

    # Concatenate slim header + all alignments, pipe into bam + sort and index
    ( cat ${prefix}_header.sam ${prefix}_alignments.sam ) \\
    | samtools view -@ ${task.cpus-1} -b - \\
    | samtools sort -@ ${task.cpus-1} -o ${output_file}
    samtools index ${output_file}

    # flagstat final file
    samtools flagstat -@ ${task.cpus-1} ${output_file} -O tsv > ${output_flagstat}

    # remove intermediates
    rm -f ${prefix}_header.sam ${prefix}_alignments.sam ${prefix}_mapped_refs.txt
    """
}

//     samtools \\
//         view \\
//         --threads ${task.cpus-1} \\
//         $args \\
//         -q $quality \\
//         -h \\
//         $input \\
//         | awk 'BEGIN{OFS="\t"} /^@/ {print; next} { split(\$3, a, "|"); if (a[2] != 9606) print }' \\
//         samtools view -b -o ${output_file} -
    
//     cat <<-END_VERSIONS > versions_samtools_view.yml
//     "${task.process}":
//         samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
//     END_VERSIONS
//     """
// }
