process KRAKEN2_PARSE {
    tag "$meta.id"

    input:
    tuple val(meta) , path(query)
    val threshold

    output:
    tuple val(meta), path("*_kraken2_filt_*"), emit: screen

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threshold = threshold
    def outfile = "${prefix}_kraken2_filt_${threshold}.txt"
    """
    krakenKmerPrcSum.sh $query $threshold > $outfile
    """
}