process KRAKEN2_PARSE {
    tag "$meta.id"

    input:
    tuple val(meta) , path(query)
    val threshold
    val kraken2_db_name

    output:
    tuple val(meta), path("*_reassign_*.txt"), emit: screen
    tuple val(meta), path("*_nReadByTaxID_preReassign*.tsv"), emit: stat_pre
    tuple val(meta), path("*_nReadByTaxID_postReassign*.tsv"), emit: stat_post

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${kraken2_db_name}"
    def script_dir = file("bin")
    def kraken_summary_script_post = script_dir.resolve("krakenFiltReadsByTaxID.py").toString()
    def kraken_summary_script_pre = script_dir.resolve("krakenReadsByTaxID.py").toString()
    def threshold = threshold
    def outfile1 = "${prefix}_reassign_${threshold}.txt"
    def outfile2 = "${prefix}_nReadByTaxID_postReassign_${threshold}.tsv"
    def outfile3 = "${prefix}_nReadByTaxID_preReassign_${threshold}.tsv"
     
    """
    krakenKmerPrcSum.sh $query $threshold > $outfile1
    python3 $kraken_summary_script_post $outfile1 $outfile2
    python3 $kraken_summary_script_pre $query $outfile3
    """
}