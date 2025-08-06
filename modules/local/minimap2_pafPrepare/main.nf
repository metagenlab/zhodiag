process PAF_PREPARE {
    tag { sample }
    input:
    tuple val(sample), val(group), path(paf_path)

    output:
    tuple val(sample), path("${sample}_annotated.paf.tsv"), emit: paf

    script:
    """
    awk 'BEGIN { OFS="\\t" }
    {
        if (NF >= 12) {
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,"${sample}","${group}"
        }
    }' ${paf_path} > ${sample}_annotated.paf.tsv
    """
 }
