process EXTRACT_CLASSIFIED {

    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::zcat'
    container "quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"

    input:
    tuple val(meta), path(classified_fastqs)

    output:
    tuple val(meta), path("*nonhuman{.,_}*"), emit: extracted_kraken2_reads

    script:
    def prefix = meta.id
    def is_single_end = meta.single_end
    def r1 = classified_fastqs[0]
    def r2 = is_single_end ? '' : classified_fastqs[1]

    """
    # ---------------- SINGLE-END ----------------
    ${is_single_end ? """
    zcat ${r1} \
      | awk 'NR%4==1 {keep=(\$0 !~ /kraken:taxid\\|9606/)} keep' \
      | gzip -c > ${prefix}_nonhuman.fastq.gz
    """ : ""}

    # ---------------- PAIRED-END ----------------
    ${!is_single_end ? """
    # R1
    zcat ${r1} \
      | awk 'NR%4==1 {keep=(\$0 !~ /kraken:taxid\\|9606/)} keep' \
      | gzip -c > ${prefix}_R1_nonhuman.fastq.gz

    # R2
    zcat ${r2} \
      | awk 'NR%4==1 {keep=(\$0 !~ /kraken:taxid\\|9606/)} keep' \
      | gzip -c > ${prefix}_R2_nonhuman.fastq.gz
    """ : ""}
    """
}
