process KRAKEN2_PARSE {
    tag "$meta.id"

    input:
    tuple val(meta) , path(query)
    val threshold

    output:
    tuple val(meta), path("*_kraken2_filt_*.txt"), emit: screen
    tuple val(meta), path("*_kraken2_filt_*_nonHuman.txt"), emit: filter
    tuple val(meta), path("*_kraken2_kmerPercByTaxID_postFiltering*.txt"), emit: report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def script_dir = file("bin")
    def kraken_summary_script = script_dir.resolve("kraken2_summary.py").toString()
    def threshold = threshold
    def outfile1 = "${prefix}_kraken2_filt_${threshold}.txt"
    def outfile2 = "${prefix}_kraken2_filt_${threshold}_nonHuman.txt"
    def outfile3 = "${prefix}_kraken2_kmerPercByTaxID_postFiltering_${threshold}.txt"
    """
    krakenKmerPrcSum.sh $query $threshold > $outfile1
    extract_nonhuman.sh $outfile1 $query > $outfile2
    python3 $kraken_summary_script $outfile2 > $outfile3
    """
}