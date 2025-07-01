process KRAKEN2_PARSE {
    tag "$meta.id"

    input:
    tuple val(meta) , path(query)
    val threshold
    val kraken2_db_name

    output:
    tuple val(meta), path("*_reassign_*.txt"), emit: screen

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${kraken2_db_name}"
    def threshold = threshold
    def outfile1 = "${prefix}_reassign_${threshold}.txt"

    """
    krakenKmerPrcSum.sh $query $threshold > $outfile1
    """
}