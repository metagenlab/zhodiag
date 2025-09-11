process SUMMARY_MAP_CANDIDATES {
    tag "$meta.id"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"


    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_summary_map_nonHuman.tsv"), emit: nonHuman

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_summary_map.tsv"
    def output_nonHuman = "${prefix}_summary_map_nonHuman.tsv"

    """
    samtools view -o ${prefix}_alignments.sam $input
    cut -f 3 ${prefix}_alignments.sam | sort | uniq -c > ${output}
    awk '{ split(\$2, a, "|"); if (a[2] != "9606") print \$1, \$2 }' OFS='\t' "${output}" > "${output_nonHuman}"
    """
}
